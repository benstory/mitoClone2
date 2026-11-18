/**************************************************************************************
 * bam2R10x.cpp An interface for R to count nucleotides in a multiplexed.bam file
 * or .cram alignment
 * The majoriy of the code was originally written for the deepSNV by Moritz Gerstung
 * Copyright 2021 story.benjamin@gmail.com
 **************************************************************************************/
/**********************************************************************
 * ORIGINAL LICENSE:
 * bamcram2R.cpp An interface for R to count nucleotides in a .bam
 * or .cram alignment
 * Copyright (C) 2015 drjsanger@github
 ***********************************************************************/
  
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
#include <string.h>
#include "htslib/sam.h"
#include "htslib/khash.h"
#include <map>
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
  
  using namespace std;
  KHASH_MAP_INIT_STR(strh,uint8_t) //readname -> readbase map used to prevent overlapping reads calls happening twice
    
    //const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";
    
    typedef struct {
      int beg, end, q, s, head_clip;
      int i;
      int* counts;
      std::map<char, int> nt_idx;
      htsFile *in;
    } nttable_t;
  
  char NUCLEOTIDES[] = {'A','T','C','G','*','N','+','-','^','$','Q'};
  int N = 11;
  
  
  std::map<std::string, int> count_valid_barcodes(const char* bamfile, const char* ref, int beg, int end, int min_count) {
    std::map<std::string, int> barcode_counts;
    htsFile *in = hts_open(bamfile, "r");
    if (!in) {
      Rf_error("Failed to open BAM file for barcode counting");
      return barcode_counts;
    }
    
    bam_hdr_t *head = sam_hdr_read(in);
    if (!head) {
      hts_close(in);
      Rf_error("Failed to read BAM header for barcode counting");
      return barcode_counts;
    }
    
    bam1_t *b = bam_init1();
    int total_reads = 0;
    int reads_with_barcodes = 0;
    
    if (strcmp(ref, "") == 0) {
      while(sam_read1(in, head, b) >= 0) {
        total_reads++;
        uint8_t *aux = bam_aux_get(b, "CB");
        if (aux) {
          char* rg = bam_aux2Z(aux);
          if (rg) {
            barcode_counts[std::string(rg)]++;
            reads_with_barcodes++;
          }
        }
      }
    } else {
      hts_idx_t *idx = sam_index_load(in, bamfile);
      if (!idx) {
        bam_destroy1(b);
        bam_hdr_destroy(head);
        hts_close(in);
        Rf_error("Failed to load BAM index for barcode counting");
        return barcode_counts;
      }
      
      //int tid = bam_name2id(head, ref);
      // std::string instead of sprintf: no fixed size buffer (R CMD check warns about sprintf)
      std::string region = std::string(ref) + ":" + std::to_string(beg) + "-" + std::to_string(end);
      hts_itr_t *iter = sam_itr_querys(idx, head, region.c_str());
      
      int result;
      while ((result = sam_itr_next(in, iter, b)) >= 0) {
        total_reads++;
        uint8_t *aux = bam_aux_get(b, "CB");
        if (aux) {
          char* rg = bam_aux2Z(aux);
          if (rg) {
            barcode_counts[std::string(rg)]++;
            reads_with_barcodes++;
          }
        }
      }

      sam_itr_destroy(iter);
      hts_idx_destroy(idx);
    }

    // Filter map to only keep barcodes with sufficient counts
    std::map<std::string, int> valid_barcodes;
    int barcodes_over_threshold = 0;
    for (const auto& pair : barcode_counts) {
      if (pair.second >= min_count) {
        valid_barcodes[pair.first] = pair.second;
        barcodes_over_threshold++;
      }
    }
    
    // Print diagnostic information
    Rprintf("FIRST PASS");
    Rprintf("Total reads processed: %d\n", total_reads);
    Rprintf("Reads with barcodes: %d (%.2f%%)\n", 
            reads_with_barcodes, 
            (reads_with_barcodes * 100.0) / total_reads);
    Rprintf("Unique barcodes found: %d\n", (int)barcode_counts.size());
    Rprintf("Barcodes with >= %d reads: %d\n", min_count, barcodes_over_threshold);
    
    bam_destroy1(b);
    bam_hdr_destroy(head);
    hts_close(in);
    
    return valid_barcodes;
  }
  
  // max.mismatches filter, checked for every read: reads without an NM tag are kept
  // (previously only the first read of the file was checked for an NM tag, which crashed R
  // when a later read had none, or switched the filter off when the first read had none)
  static bool nm_ok(const bam1_t *b, int maxmismatches) {
    if (maxmismatches == -1) return true;
    uint8_t *nm = bam_aux_get(b, "NM");
    return !nm || bam_aux2i(nm) <= maxmismatches;
  }
  

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
  
  extern "C" {
    void add_new_tag(std::map<std::string, std::vector<int> > *dict, const char* new_key, int beg, int end)
    {
      // //cout<< "Added tag: " << new_key << endl;
      std::vector<int> counts((end-beg+1)*11*2,0);
      dict->insert({new_key, counts});
    }
    void bam2R_pileup_function(const bam_pileup1_t *pl, int pos, int n_plp, nttable_t& nttable,std::map<std::string,std::vector<int> >& dic,const std::map<std::string, int>& valid_barcodes)
    {
      int i, s;
      // needs to be moved
      int missingBC = 0;
      int len = nttable.end - nttable.beg;
      std::map<char, int> nt_freq;
      khash_t(strh) *h;
      khiter_t k;
      h = kh_init(strh);
      
      //if current position is greater than the start and less than the end position
      if ((int)pos >= nttable.beg && (int)pos < nttable.end)
      {
        //set relative counts pointer to distance to start pos
        //int* counts = nttable.counts + (int)pos - nttable.beg ;
        int count_offset = (int)pos - nttable.beg ;
        //iterate through i starting at 0 to n_plp	
        for (i=0; i<n_plp; i++)
        {
          const bam_pileup1_t *p = pl + i;
          // check if the const char* is reverse
          s = bam_is_rev(p->b) * len * N;
          int absent;
          k = kh_put(strh, h, bam_get_qname(p->b), &absent);
          uint8_t cbase = bam_seqi(bam_get_seq(p->b),p->qpos);	    
          uint8_t pre_b;
          
          //cleaned barcode reading
          //const char *rg;
          //rg = (char *)bam_aux_get(p->b, "CB");	 
          uint8_t *aux = bam_aux_get(p->b, "CB");
          if (!aux) {
            missingBC++;
            continue;
          }
          char* rg = bam_aux2Z(aux);
          if (!rg) {
            Rf_warning("Barcode is missing in some reads!");
            missingBC++;
          }else{
            if (!valid_barcodes.count(std::string(rg))) {
              continue; // Skip this read if barcode not in valid set
            }
            if (dic.count(rg) == 0){
              add_new_tag(&dic, rg, nttable.beg, nttable.end);
            }	   
            if(!absent){ //Read already processed to get base processed (we only increment if base is different between overlapping read pairs)
              k = kh_get(strh, h, bam_get_qname(p->b));
              pre_b = kh_val(h,k);
            }else{
              //Add the value to the hash
              kh_value(h, k) = cbase;
            }
            {      
              if(!absent && pre_b == cbase) continue;	     	     
              if (p->is_tail){
                //counts[s + len * nttable.nt_idx['$']]++;
                dic[rg][count_offset + s + len * nttable.nt_idx['$']]++;
              }
              else if (p->is_head){
                //counts[s + len * nttable.nt_idx['^']]++;
                dic[rg][count_offset + s + len * nttable.nt_idx['^']]++;
              }
              // distance to the 5' end (head) of the read: for reverse reads that is the END of the
              // stored sequence (previously the wrong end of reverse reads was clipped, one base too few)
              int from_head = bam_is_rev(p->b) ? p->b->core.l_qseq - 1 - p->qpos : p->qpos;
              if(from_head < nttable.head_clip){
                //counts[s + len * nttable.nt_idx['N']]++;
                //ADD N's
                dic[rg][count_offset + s + len * nttable.nt_idx['N']]++;
              }else{
                if (!p->is_del) {
                  int  c = seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)];
                  if( bam_get_qual(p->b)[p->qpos] > nttable.q){
                    //counts[s + len * nttable.nt_idx[char(c)]]++;
                    dic[rg][count_offset + s + len * nttable.nt_idx[char(c)]]++;
                  }
                  else{
                    //counts[s + len * nttable.nt_idx['N']]++;
                    dic[rg][count_offset + s + len * nttable.nt_idx['N']]++;
                  }
                  if (p->indel > 0){
                    //counts[s + len * nttable.nt_idx['+']]++;
                    dic[rg][count_offset + s + len * nttable.nt_idx['+']]++;
                  }
                  else if (p->indel < 0){
                    //counts[s + len * nttable.nt_idx['-']]++;
                    dic[rg][count_offset + s + len * nttable.nt_idx['-']]++;
                  }
                } else{
                  //counts[s + len * nttable.nt_idx['*']]++;
                  dic[rg][count_offset + s + len * nttable.nt_idx['*']]++;
                }
                //counts[s + len * nttable.nt_idx['Q']] += p->b->core.qual;
                dic[rg][count_offset + s + len * nttable.nt_idx['Q']] += p->b->core.qual;
              }
            }
          }
        }
        nttable.i++;
      }
      kh_destroy(strh, h);
    }  
    SEXP bam2R_10x(SEXP bamfiletest, SEXP reftest, SEXP begtest, SEXP endtest, int countstest[], SEXP qtest, SEXP mqtest, SEXP stest, SEXP head_cliptest, SEXP maxdepthtest, SEXP verbosetest, SEXP masktest, SEXP keepflagtest, SEXP maxmismatchestest, SEXP min_reads_per_barcodetest)
    {
      const char* bamfile = CHAR(STRING_ELT(bamfiletest,0));
      const char* ref = CHAR(STRING_ELT(reftest,0));
      int* beg = INTEGER(begtest);
      int* end = INTEGER(endtest);
      //int counts[(*end-*beg+1)*11*2];
      int* q = INTEGER(qtest);
      int* mq = INTEGER(mqtest);
      int* s = INTEGER(stest);
      int* head_clip = INTEGER(head_cliptest);
      int* maxdepth = INTEGER(maxdepthtest);
      int* verbose = INTEGER(verbosetest);
      int* mask = INTEGER(masktest);
      int* keepflag = INTEGER(keepflagtest);
      int* maxmismatches = INTEGER(maxmismatchestest);
      int* min_reads_per_barcode = INTEGER(min_reads_per_barcodetest);  
      
      bam_plp_t buf = NULL;
      bam1_t *b = NULL;
      //hts_itr_t *iter = NULL;
      bam_hdr_t *head = NULL;
      
      //int c = 0;
      nttable_t nttable;
      nttable.q = *q; //Base quality cutoff
      nttable.s = *s; //Strand (2=both)
      nttable.head_clip = *head_clip;
      nttable.i = 0;
      
      //valid barcodes identifier
      auto valid_barcodes = count_valid_barcodes(bamfile, ref, *beg, *end, *min_reads_per_barcode);
      
      std::map<std::string,std::vector<int> > dic;
      //nttable.counts = counts;
      int i;
      for (i=0; i<N; i++)
        nttable.nt_idx[NUCLEOTIDES[i]] = i;
      nttable.beg = *beg -1;
      nttable.end = *end;
      nttable.in = hts_open(bamfile, "r");
      if (nttable.in == 0) {
        Rf_error("Fail to open input BAM/CRAM file %s\n", bamfile);
        return Rf_ScalarLogical(FALSE);//return 1;
      }
      buf = bam_plp_init(0,(void *)&nttable); // initialize pileup
      bam_plp_set_maxcnt(buf,*maxdepth);
      b = bam_init1();
      //get header
      head = sam_hdr_read(nttable.in);
      //int mask = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY;
      int tid, pos, n_plp = -1;
      const bam_pileup1_t *pl;
      //get first read to check for NM tag
      int ret;
      ret = sam_read1(nttable.in, head, b);
      hts_close(nttable.in);
      nttable.in = hts_open(bamfile, "r");
      uint8_t *paux = bam_aux_get(b, "NM");
      if ( ! paux && *maxmismatches != -1 ) {
        Rf_warning("The first read has no NM tag: reads without an NM tag are not filtered by max.mismatches.\n");
        //return Rf_ScalarLogical(FALSE);
      }
      if (strcmp(ref, "") == 0) { // if a region is not specified
        //Replicate sampileup functionality (uses above mask without supplementary)
        while((ret = sam_read1(nttable.in, head, b)) >= 0){
          if ((b->core.flag & *mask)==0 && b->core.qual >= *mq && (b->core.flag & *keepflag)==*keepflag && nm_ok(b, *maxmismatches)){	 
            bam_plp_push(buf, b);
          };
          while ( (pl=bam_plp_next(buf, &tid, &pos, &n_plp)) != 0) {
            bam2R_pileup_function(pl,pos,n_plp,nttable,dic,valid_barcodes);
          }
        }
      }
      else {
        int tid;
        hts_idx_t *idx;
        idx = sam_index_load(nttable.in,bamfile); // load BAM index
        if (idx == 0) {
          Rf_error("BAM/CRAM index file is not available.\n");
          return Rf_ScalarLogical(FALSE);//return 1;
        }
        tid = bam_name2id(head, ref);
        if (tid < 0) {
          char msg[1024];
          invalid_sequence_message(msg, sizeof(msg), ref, head);
          Rf_error("%s\n", msg);
          return Rf_ScalarLogical(FALSE);//return 1;
        }
        // std::string instead of sprintf: no fixed size buffer (R CMD check warns about sprintf)
        std::string region = std::string(ref) + ":" + std::to_string(nttable.beg) + "-" + std::to_string(nttable.end);
        if(*verbose){
          Rprintf("Reading %s at coordinates %s.\n", bamfile, region.c_str());
        }
        //Implement a fetch style iterator
        hts_itr_t *iter = sam_itr_querys(idx, head, region.c_str());
        int result;
        while ((result = sam_itr_next(nttable.in, iter, b)) >= 0) {
          if ((b->core.flag & *mask)==0 && b->core.qual >= *mq && (b->core.flag & *keepflag)==*keepflag && nm_ok(b, *maxmismatches)){
            bam_plp_push(buf, b);
          };
          while ( (pl=bam_plp_next(buf, &tid, &pos, &n_plp)) != 0) {
            bam2R_pileup_function(pl,pos,n_plp,nttable,dic,valid_barcodes);
          }
        }
        if(result < -1){
          Rf_error("Error code (%d) encountered reading sam iterator.\n", result);
          return Rf_ScalarLogical(FALSE);//return 1;
        }
        sam_itr_destroy(iter);
        hts_idx_destroy(idx);
      }

      bam_plp_push(buf,0); // finalize pileup
      
      while ( (pl=bam_plp_next(buf, &tid, &pos, &n_plp)) != 0) {
        bam2R_pileup_function(pl,pos,n_plp,nttable, dic,valid_barcodes);
      }
      bam_destroy1(b);
      bam_hdr_destroy(head);
      bam_plp_destroy(buf);
      hts_close(nttable.in);  
      //return 0;
      // how many barcodes exist
      int listSize = dic.size();
      if(*verbose){
        Rprintf("Detected %d unique barcodes.\n", listSize);
      }
      // initialize iterator
      int key_index = 0;
      // protect R output - names of list and values in list
      SEXP labels = PROTECT(Rf_allocVector(STRSXP, listSize));
      SEXP vec = PROTECT(Rf_allocVector(VECSXP, listSize));
      // how long is each vector element in the list
      int vecSize = (*end-*beg+1)*11*2;
      // pop whole dictionary and create names
      for (auto const& pos : dic)
      {
        // pull out the key
        std::string keyinit = pos.first;
        const char* key = keyinit.c_str(); 
        // assign the key to an index that contains a numeric vector of length vecSize
        SET_VECTOR_ELT(vec, key_index, Rf_allocVector(INTSXP, vecSize));
        // store the names of the keys at the same index in a string vectors
        SET_STRING_ELT(labels, key_index, Rf_mkChar(key));
        // iterate over counts for this label and set elements
        for (int i = 0; i < vecSize; i++){
          INTEGER(VECTOR_ELT(vec, key_index))[i] = pos.second[i];
        }
        key_index++;
      }
      //assign names to list
      Rf_setAttrib(vec, R_NamesSymbol, labels);
      //cleanup and return
      UNPROTECT(2);
      return vec;
    }
    
    R_CMethodDef cMethods[] = {
      {"bam2R_10x", (DL_FUNC) &bam2R_10x, 15}
    };
    
    void R_init_bam2R_10x(DllInfo *info) {
      R_registerRoutines(info, cMethods, NULL, NULL, NULL);
    }
  } // extern "C"
