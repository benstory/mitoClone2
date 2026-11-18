/**************************************************************************************
 * bam2R_v2.cpp  Version 2 of bam2R_10x.cpp: count nucleotides per cell barcode in a
 * multiplexed .bam or .cram alignment
 * The majoriy of the code was originally written for the deepSNV by Moritz Gerstung
 * Copyright 2021 story.benjamin@gmail.com
 *
 * This file is INDEPENDENT of bam2R_10x.cpp (version 1). Both can be compiled into
 * the same shared library:
 * - every helper lives in an anonymous namespace (no symbol clashes with version 1,
 *   e.g. N, nttable_t or bam2R_pileup_function exist in both files)
 * - the only exported symbols are the two .Call entry points
 *     bam2R_v2_barcodes  (first pass, reads per barcode)
 *     bam2R_v2_count     (sparse nucleotide counts)
 *     bam2R_v2_fasta     (reference sequence of regions from a FASTA file)
 *   they are called by bam2R_v2() in R/bam_functions_v2.R
 *
 * SPARSE VERSION:
 * - any number of regions (any chromosome) are read in ONE pass over the file
 * - counts are never stored as a dense (positions x 22) vector per barcode, only
 *   the (site, barcode, column) combinations that were actually observed are kept
 * - the barcode of a read is looked up once per read, not once per read per position
 * - both strands are reported separately (needed for the strand bias)
 * - optional UMI mode: a UMI of a cell is counted once per site/nucleotide/strand
 *   instead of counting every read (PCR duplicates)
 **************************************************************************************/
/**********************************************************************
 * ORIGINAL LICENSE:
 * bamcram2R.cpp An interface for R to count nucleotides in a .bam
 * or .cram alignment
 * Copyright (C) 2015 drjsanger@github
 ***********************************************************************/

#include <vector>
#include <string>
#include <algorithm>
#include <exception>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/khash.h"
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

namespace {

  KHASH_MAP_INIT_STR(v2_strh,uint8_t) //readname -> readbase map used to prevent overlapping reads calls happening twice
  KHASH_MAP_INIT_STR(v2_bch,int)      //barcode -> barcode index (whitelist) or barcode -> read count (first pass)
  KHASH_SET_INIT_INT64(v2_umih)       //(UMI, barcode, column) already counted at the current position

  // columns per strand, forward = 0..7 and reverse = 8..15
  // R names:  A  T  C  G  -  N  INS DEL   a  t  c  g  _  n  ins del
  const int N = 8;
  enum {IDX_GAP = 4, IDX_N = 5, IDX_INS = 6, IDX_DEL = 7};
  // htslib 4-bit base (=ACMGRSVTWYHKDBN) -> column, everything that is not A/C/G/T is an N
  const int NT16_IDX[16] = {5,0,2,5, 3,5,5,5, 1,5,5,5, 5,5,5,5};

  // a region of interest, 0-based half open
  // offset = number of requested positions in front of this region (provided by R)
  // every requested position therefore has its own site number: offset + distance to beg + 1
  typedef struct {
    int tid, beg, end, offset;
  } region_t;

  bool region_cmp(const region_t& a, const region_t& b) {
    return a.tid < b.tid || (a.tid == b.tid && a.beg < b.beg);
  }

  // everything needed to walk a set of regions in a bam file
  typedef struct {
    htsFile *in;
    bam_hdr_t *head;
    hts_idx_t *idx;
    hts_itr_t *iter;
  } bamreader_t;

  typedef struct {
    int q, mq, head_clip, mask, keepflag, maxmismatches;
    char tag[3];
    char umitag[3];                // UMI tag, empty = count reads instead of UMIs
    uint32_t cur_umi;              // UMI (hashed) of the read that is being added to the pileup
    khash_t(v2_umih) *umis;        // (UMI, barcode, column) seen at the current position
    int status;                    // negative if reading the file failed
    int cur_bc;                    // barcode index of the read that is being added to the pileup
    size_t cur_reg;                // region the pileup is currently in
    bamreader_t reader;
    std::vector<region_t> regions; // sorted by tid, beg
    khash_t(v2_bch) *barcodes;     // valid barcodes -> index
    khash_t(v2_strh) *seen;        // read names seen at the current position
    std::vector<int> counts;       // counts of the current position, one per (barcode, column)
    std::vector<int> keys;         // (barcode, column) that have a count at the current position
    // OUTPUT: the covered sites ...
    std::vector<int> sites;
    // ... and for each of the 16 columns one entry per observed (site, barcode)
    std::vector<int> o_site[16], o_bc[16], o_count[16];
  } nttable_t;

  char errbuf[1024];

  // "Invalid sequence X" plus the sequence names of the file, so that a chrM / MT mixup is
  // obvious. Written into buf (size bytes).
  static void invalid_sequence_message(char* buf, size_t size, const char* ref, const bam_hdr_t* head) {
    size_t n = (size_t) snprintf(buf, size, "Invalid sequence %s. ", ref);
    const char* mito = NULL;
    for (int t = 0; t < head->n_targets; t++) {
      const char* name = head->target_name[t];
      if (!strcmp(name, "chrM") || !strcmp(name, "MT") || !strcmp(name, "M") || !strcmp(name, "chrMT")) mito = name;
    }
    if (mito && n < size)
      n += (size_t) snprintf(buf + n, size - n, "The mitochondrial genome is called %s in this file (e.g. sites = \"%s:1-16569\"). ", mito, mito);
    if (n < size) n += (size_t) snprintf(buf + n, size - n, "The file contains: ");
    for (int t = 0; t < head->n_targets && t < 10 && n < size; t++)
      n += (size_t) snprintf(buf + n, size - n, "%s%s", t ? ", " : "", head->target_name[t]);
    if (head->n_targets > 10 && n < size)
      snprintf(buf + n, size - n, ", ... (%d sequences)", head->n_targets);
  }

  void close_bamreader(bamreader_t& r) {
    if (r.iter) hts_itr_destroy(r.iter);
    if (r.idx) hts_idx_destroy(r.idx);
    if (r.head) bam_hdr_destroy(r.head);
    if (r.in) hts_close(r.in);
    r.iter = NULL; r.idx = NULL; r.head = NULL; r.in = NULL;
  }

  // open the file + index and set up ONE iterator that covers all regions
  // returns 0 if everything worked, otherwise errbuf contains the reason
  int open_bamreader(bamreader_t& r, const char* bamfile, SEXP chr, SEXP start, SEXP end, SEXP offset, std::vector<region_t>* regions)
  {
    r.in = NULL; r.head = NULL; r.idx = NULL; r.iter = NULL;
    int nreg = Rf_length(chr);
    r.in = hts_open(bamfile, "r");
    if (r.in == 0) {
      snprintf(errbuf, sizeof(errbuf), "Fail to open input BAM/CRAM file %s", bamfile);
      return 1;
    }
    r.head = sam_hdr_read(r.in);
    if (r.head == 0) {
      snprintf(errbuf, sizeof(errbuf), "Fail to read the header of %s", bamfile);
      return 1;
    }
    r.idx = sam_index_load(r.in, bamfile); // load BAM index
    if (r.idx == 0) {
      snprintf(errbuf, sizeof(errbuf), "BAM/CRAM index file is not available.");
      return 1;
    }
    // region strings for htslib, {} protects contig names that contain a colon
    std::vector<std::string> regstr(nreg);
    std::vector<char*> regarray(nreg);
    for (int i = 0; i < nreg; i++) {
      const char* ref = CHAR(STRING_ELT(chr, i));
      int tid = bam_name2id(r.head, ref);
      if (tid < 0) {
        invalid_sequence_message(errbuf, sizeof(errbuf), ref, r.head);
        return 1;
      }
      regstr[i] = "{" + std::string(ref) + "}:" + std::to_string(INTEGER(start)[i]) + "-" + std::to_string(INTEGER(end)[i]);
      regarray[i] = &regstr[i][0];
      if (regions) {
        region_t reg = {tid, INTEGER(start)[i] - 1, INTEGER(end)[i], INTEGER(offset)[i]};
        regions->push_back(reg);
      }
    }
    if (regions) std::sort(regions->begin(), regions->end(), region_cmp);
    //Implement a fetch style iterator over all regions at once
    r.iter = sam_itr_regarray(r.idx, r.head, regarray.data(), nreg);
    if (r.iter == 0) {
      snprintf(errbuf, sizeof(errbuf), "Fail to set up the iterator for the provided regions.");
      return 1;
    }
    return 0;
  }

  // UMI sequence -> 32 bit number. UMIs of up to 15 A/C/G/T (10x: 10 or 12 bases) are packed
  // with 2 bits per base behind a leading 1, which is lossless: two different UMIs never get
  // the same number. Everything else is hashed into the upper half (highest bit set), so the
  // two kinds can not collide either.
  uint32_t encode_umi(const char* umi)
  {
    uint32_t v = 1;
    int n = 0;
    for (const char* c = umi; *c; c++, n++) {
      uint32_t code;
      switch (*c) {
      case 'A': code = 0; break;
      case 'C': code = 1; break;
      case 'G': code = 2; break;
      case 'T': code = 3; break;
      default: return (uint32_t) kh_str_hash_func(umi) | 0x80000000u;
      }
      if (n == 15) return (uint32_t) kh_str_hash_func(umi) | 0x80000000u;
      v = (v << 2) | code;
    }
    return v;
  }

  // hands the next USABLE read to the pileup, all read level filters happen here
  // so reads of unwanted barcodes never enter the pileup at all
  int read_function(void *data, bam1_t *b)
  {
    nttable_t *nttable = (nttable_t*) data;
    int result;
    while ((result = sam_itr_next(nttable->reader.in, nttable->reader.iter, b)) >= 0) {
      if ((b->core.flag & nttable->mask) != 0 || b->core.qual < nttable->mq || (b->core.flag & nttable->keepflag) != nttable->keepflag) continue;
      if (nttable->maxmismatches != -1) {
        // reads without an NM tag are kept
        uint8_t *paux = bam_aux_get(b, "NM");
        if (paux && bam_aux2i(paux) > nttable->maxmismatches) continue;
      }
      // the barcode index of the read. An empty tag means that the whole file is one cell
      // (one BAM file per cell, baseCountsFromBamList_v2): every read gets index 0
      int bc = 0;
      if (nttable->tag[0] != '\0') {
        uint8_t *aux = bam_aux_get(b, nttable->tag);
        if (!aux) continue;
        const char* rg = bam_aux2Z(aux);
        if (!rg) continue;
        khiter_t k = kh_get(v2_bch, nttable->barcodes, rg);
        if (k == kh_end(nttable->barcodes)) continue; // Skip this read if barcode not in valid set
        bc = kh_value(nttable->barcodes, k);
      }
      if (nttable->umitag[0] != '\0') {
        // reads without a UMI can not be collapsed and are skipped
        uint8_t *umiaux = bam_aux_get(b, nttable->umitag);
        if (!umiaux) continue;
        const char* umi = bam_aux2Z(umiaux);
        if (!umi) continue;
        nttable->cur_umi = encode_umi(umi);
      }
      nttable->cur_bc = bc;
      break;
    }
    if (result < -1) nttable->status = result;
    return result;
  }

  // called once when a read enters the pileup: remember its barcode index (lower 32 bit)
  // and its UMI (upper 32 bit)
  int pileup_constructor(void *data, const bam1_t *b, bam_pileup_cd *cd)
  {
    // unsigned arithmetic: shifting a signed value into the sign bit is undefined behaviour
    nttable_t *nttable = (nttable_t*) data;
    cd->i = (int64_t) ((((uint64_t) nttable->cur_umi) << 32) | (uint32_t) nttable->cur_bc);
    return 0;
  }

  void pileup_function(const bam_pileup1_t *pl, int tid, int pos, int n_plp, nttable_t& nttable)
  {
    int i, s;
    //move on to the first region that does not end before the current position
    while (nttable.cur_reg < nttable.regions.size() &&
           (nttable.regions[nttable.cur_reg].tid < tid ||
            (nttable.regions[nttable.cur_reg].tid == tid && nttable.regions[nttable.cur_reg].end <= pos)))
      nttable.cur_reg++;
    if (nttable.cur_reg == nttable.regions.size()) return;
    const region_t& reg = nttable.regions[nttable.cur_reg];
    //only continue if the current position is inside of a region
    if (reg.tid != tid || pos < reg.beg) return;

    khash_t(v2_strh) *h = nttable.seen;
    khiter_t k;
    kh_clear(v2_strh, h);
    khash_t(v2_umih) *u = nttable.umis;
    bool use_umi = nttable.umitag[0] != '\0';
    if (use_umi) kh_clear(v2_umih, u);
    int umi_absent;
    int* counts = nttable.counts.data();
    std::vector<int>& keys = nttable.keys;
    // 16 counts per barcode: a count that is hit for the first time is remembered in keys
    // UMI mode: every (UMI, barcode, column) is only counted once per position
#define ADD_COUNT(col) { int key = bc + s + (col); \
      if (use_umi) kh_put(v2_umih, u, (umi << 32) | (uint32_t) key, &umi_absent); \
      if (!use_umi || umi_absent) { if (counts[key]++ == 0) keys.push_back(key); } }
    //iterate through i starting at 0 to n_plp
    for (i=0; i<n_plp; i++)
    {
      const bam_pileup1_t *p = pl + i;
      // reverse strand counts are shifted by N columns
      s = bam_is_rev(p->b) * N;
      int bc = ((int) (p->cd.i & 0xffffffff)) * 2 * N;
      uint64_t umi = ((uint64_t) p->cd.i) >> 32;
      if (p->b->core.flag & BAM_FPAIRED) {
        //Read already processed to get base processed (we only increment if base is different between overlapping read pairs)
        int absent;
        uint8_t cbase = bam_seqi(bam_get_seq(p->b),p->qpos);
        k = kh_put(v2_strh, h, bam_get_qname(p->b), &absent);
        if (!absent && kh_val(h,k) == cbase) continue;
        if (absent) kh_value(h, k) = cbase; //Add the value to the hash
      }
      // distance to the 5' end (head) of the read: for reverse reads that is the END of the
      // stored sequence. (The deepSNV code clipped the wrong end of reverse reads and one
      // base too few: qpos < head_clip || l_qseq - qpos < head_clip.)
      int from_head = bam_is_rev(p->b) ? p->b->core.l_qseq - 1 - p->qpos : p->qpos;
      if (from_head < nttable.head_clip) {
        //ADD N's
        ADD_COUNT(IDX_N);
      }else{
        if (!p->is_del) {
          if( bam_get_qual(p->b)[p->qpos] > nttable.q){
            ADD_COUNT(NT16_IDX[bam_seqi(bam_get_seq(p->b), p->qpos)]);
          }
          else{
            ADD_COUNT(IDX_N);
          }
          if (p->indel > 0){
            ADD_COUNT(IDX_INS);
          }
          else if (p->indel < 0){
            ADD_COUNT(IDX_DEL);
          }
        } else{
          ADD_COUNT(IDX_GAP);
        }
      }
    }
#undef ADD_COUNT
    if (keys.empty()) return;
    //store the non-zero counts of this site and reset them for the next position
    int site = reg.offset + (pos - reg.beg) + 1;
    nttable.sites.push_back(site);
    for (size_t j = 0; j < keys.size(); j++) {
      int col = keys[j] % (2 * N);
      nttable.o_site[col].push_back(site);
      nttable.o_bc[col].push_back(keys[j] / (2 * N) + 1);
      nttable.o_count[col].push_back(counts[keys[j]]);
      counts[keys[j]] = 0;
    }
    keys.clear();
  }

  // copy a vector to R and release its memory right away
  SEXP int_vector_to_R(std::vector<int>& v)
  {
    SEXP out = PROTECT(Rf_allocVector(INTSXP, v.size()));
    if (!v.empty()) memcpy(INTEGER(out), v.data(), v.size() * sizeof(int));
    std::vector<int>().swap(v);
    UNPROTECT(1);
    return out;
  }

  // the same for the 16 columns: list of 16 integer vectors
  SEXP int_columns_to_R(std::vector<int>* v)
  {
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 2 * N));
    for (int i = 0; i < 2 * N; i++)
      SET_VECTOR_ELT(out, i, int_vector_to_R(v[i]));
    UNPROTECT(1);
    return out;
  }

} // anonymous namespace

extern "C" {

  // FIRST PASS (optional): number of reads per barcode within the regions
  // only reads the aux tag, no pileup. Skipped completely if the user provides barcodes
  SEXP bam2R_v2_barcodes(SEXP bamfiletest, SEXP chrtest, SEXP starttest, SEXP endtest, SEXP tagtest)
  {
    const char* bamfile = CHAR(STRING_ELT(bamfiletest,0));
    const char* tag = CHAR(STRING_ELT(tagtest,0));
    bamreader_t reader;
    if (open_bamreader(reader, bamfile, chrtest, starttest, endtest, R_NilValue, NULL)) {
      close_bamreader(reader);
      Rf_error("%s", errbuf);
    }
    khash_t(v2_bch) *h = kh_init(v2_bch);
    khiter_t k;
    bam1_t *b = bam_init1();
    int result, absent;
    while ((result = sam_itr_next(reader.in, reader.iter, b)) >= 0) {
      uint8_t *aux = bam_aux_get(b, tag);
      if (!aux) continue;
      const char* rg = bam_aux2Z(aux);
      if (!rg) continue;
      k = kh_get(v2_bch, h, rg);
      if (k == kh_end(h)) {
        // the key has to outlive the read
        k = kh_put(v2_bch, h, strdup(rg), &absent);
        kh_value(h, k) = 0;
      }
      kh_value(h, k)++;
    }
    bam_destroy1(b);
    close_bamreader(reader);
    // named integer vector: barcode -> number of reads
    int listSize = kh_size(h);
    SEXP labels = PROTECT(Rf_allocVector(STRSXP, listSize));
    SEXP vec = PROTECT(Rf_allocVector(INTSXP, listSize));
    int key_index = 0;
    for (k = kh_begin(h); k != kh_end(h); ++k) {
      if (!kh_exist(h, k)) continue;
      SET_STRING_ELT(labels, key_index, Rf_mkChar(kh_key(h, k)));
      INTEGER(vec)[key_index] = kh_value(h, k);
      free((char*) kh_key(h, k));
      key_index++;
    }
    kh_destroy(v2_bch, h);
    Rf_setAttrib(vec, R_NamesSymbol, labels);
    UNPROTECT(2);
    if (result < -1) Rf_error("Error code (%d) encountered reading sam iterator.", result);
    return vec;
  }

  // SECOND PASS: sparse nucleotide counts per (site, barcode, column)
  SEXP bam2R_v2_count(SEXP bamfiletest, SEXP chrtest, SEXP starttest, SEXP endtest, SEXP offsettest, SEXP barcodestest, SEXP qtest, SEXP mqtest, SEXP head_cliptest, SEXP maxdepthtest, SEXP verbosetest, SEXP masktest, SEXP keepflagtest, SEXP maxmismatchestest, SEXP tagtest, SEXP umitagtest)
  {
    const char* bamfile = CHAR(STRING_ELT(bamfiletest,0));
    int maxdepth = INTEGER(maxdepthtest)[0];
    int verbose = INTEGER(verbosetest)[0];

    int nbc = Rf_length(barcodestest);
    nttable_t *nttable = NULL;
    bam_plp_t buf = NULL;
    // 0 = ok, 1 = the reason is in errbuf
    int failed = 0;
    // a C++ exception (std::bad_alloc when the sparse output does not fit into memory) must
    // never reach R, it would abort the whole R session: catch it and turn it into an R error
    try {
      nttable = new nttable_t;
      nttable->q = INTEGER(qtest)[0]; //Base quality cutoff
      nttable->mq = INTEGER(mqtest)[0];
      nttable->head_clip = INTEGER(head_cliptest)[0];
      nttable->mask = INTEGER(masktest)[0];
      nttable->keepflag = INTEGER(keepflagtest)[0];
      nttable->maxmismatches = INTEGER(maxmismatchestest)[0];
      strncpy(nttable->tag, CHAR(STRING_ELT(tagtest,0)), 2);
      nttable->tag[2] = '\0';
      // empty string = no UMI collapsing
      strncpy(nttable->umitag, CHAR(STRING_ELT(umitagtest,0)), 2);
      nttable->umitag[2] = '\0';
      nttable->cur_umi = 0;
      nttable->status = 0;
      nttable->cur_bc = 0;
      nttable->cur_reg = 0;
      nttable->barcodes = NULL;
      nttable->seen = NULL;
      nttable->umis = NULL;
      nttable->reader.in = NULL; nttable->reader.head = NULL; nttable->reader.idx = NULL; nttable->reader.iter = NULL;

      if (open_bamreader(nttable->reader, bamfile, chrtest, starttest, endtest, offsettest, &nttable->regions)) {
        failed = 1;
      } else {
        //valid barcodes -> index, the keys are owned by R and stay alive during this call
        int absent;
        nttable->barcodes = kh_init(v2_bch);
        for (int i = 0; i < nbc; i++) {
          khiter_t k = kh_put(v2_bch, nttable->barcodes, CHAR(STRING_ELT(barcodestest, i)), &absent);
          kh_value(nttable->barcodes, k) = i;
        }
        nttable->seen = kh_init(v2_strh);
        nttable->umis = kh_init(v2_umih);
        nttable->counts.assign((size_t) nbc * 2 * N, 0);

        // initialize pileup, it pulls its reads from read_function
        buf = bam_plp_init(read_function, (void *) nttable);
        bam_plp_constructor(buf, pileup_constructor);
        bam_plp_set_maxcnt(buf, maxdepth);
        int tid, pos, n_plp = -1;
        const bam_pileup1_t *pl;
        while ( (pl=bam_plp_auto(buf, &tid, &pos, &n_plp)) != 0) {
          pileup_function(pl, tid, pos, n_plp, *nttable);
        }
        if (nttable->status < 0) {
          snprintf(errbuf, sizeof(errbuf), "Error code (%d) encountered reading sam iterator.", nttable->status);
          failed = 1;
        }
      }
    } catch (std::exception& e) {
      snprintf(errbuf, sizeof(errbuf), "bam2R_v2 ran out of memory (%s). Try fewer barcodes, fewer sites or more ncores.", e.what());
      failed = 1;
    }
    // free everything that is not needed for the output (all of it if something failed)
    if (buf) bam_plp_destroy(buf);
    if (nttable) {
      if (nttable->seen) kh_destroy(v2_strh, nttable->seen);
      if (nttable->umis) kh_destroy(v2_umih, nttable->umis);
      if (nttable->barcodes) kh_destroy(v2_bch, nttable->barcodes);
      close_bamreader(nttable->reader);
      std::vector<int>().swap(nttable->counts);
    }
    if (failed) {
      delete nttable;
      Rf_error("%s", errbuf);
    }
    if(verbose){
      Rprintf("Detected %d covered sites.\n", (int)nttable->sites.size());
    }
    // hand the sparse counts over to R
    const char *names[] = {"sites", "site", "barcode", "count", ""};
    SEXP vec = PROTECT(Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(vec, 0, int_vector_to_R(nttable->sites));
    SET_VECTOR_ELT(vec, 1, int_columns_to_R(nttable->o_site));
    SET_VECTOR_ELT(vec, 2, int_columns_to_R(nttable->o_bc));
    SET_VECTOR_ELT(vec, 3, int_columns_to_R(nttable->o_count));
    delete nttable;
    //cleanup and return
    UNPROTECT(1);
    return vec;
  }

  // Reference sequence of regions (1-based, inclusive) from a FASTA file. The .fai index
  // is created next to the file if it does not exist. NA for sequences that are not in
  // the file, used by mutationCallsFromCohort(all.chromosomes = TRUE, reference = <file>)
  SEXP bam2R_v2_fasta(SEXP fastatest, SEXP chrtest, SEXP starttest, SEXP endtest)
  {
    const char* fasta = CHAR(STRING_ELT(fastatest, 0));
    faidx_t *fai = fai_load(fasta);
    if (fai == NULL) {
      Rf_error("Fail to open or index the FASTA file %s (the folder has to be writable to create the .fai index, or run 'samtools faidx' first).", fasta);
    }
    int n = Rf_length(chrtest);
    SEXP out = PROTECT(Rf_allocVector(STRSXP, n));
    for (int i = 0; i < n; i++) {
      const char* chr = CHAR(STRING_ELT(chrtest, i));
      SET_STRING_ELT(out, i, NA_STRING);
      if (!faidx_has_seq(fai, chr)) continue;
      // a region that starts behind the end of the sequence gets an empty string (the
      // sites have no reference base). htslib would clamp it and return the last base.
      if (INTEGER(starttest)[i] < 1 || (hts_pos_t) INTEGER(starttest)[i] > faidx_seq_len64(fai, chr)) {
        SET_STRING_ELT(out, i, R_BlankString);
        continue;
      }
      int len = 0;
      char *seq = faidx_fetch_seq(fai, chr, INTEGER(starttest)[i] - 1, INTEGER(endtest)[i] - 1, &len);
      if (seq == NULL) continue;
      if (len >= 0) SET_STRING_ELT(out, i, Rf_mkCharLen(seq, len));
      free(seq);
    }
    fai_destroy(fai);
    UNPROTECT(1);
    return out;
  }

  // No R_init_* here on purpose: the package init function (if any) lives in one place.
  // If the package registers its routines and turns dynamic lookup off, add:
  //   {"bam2R_v2_count",    (DL_FUNC) &bam2R_v2_count,    16},
  //   {"bam2R_v2_barcodes", (DL_FUNC) &bam2R_v2_barcodes,  5},
  //   {"bam2R_v2_fasta",    (DL_FUNC) &bam2R_v2_fasta,     4},

} // extern "C"
