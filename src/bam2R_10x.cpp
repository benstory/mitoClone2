/**************************************************************************************
 * bam2R10x.cpp An interface for R to count nucleotides in a multiplexed.bam file
 * or .cram alignment
 * The majoriy of the code was originally written for the deepSNV by Moritz Gerstung
 * Copyright 2021 story.benjamin@gmail.com
 /**************************************************************************************
#ORIGINAL LICENSE:
/**********************************************************************
 * bamcram2R.cpp An interface for R to count nucleotides in a .bam
 * or .cram alignment
 * Copyright (C) 2015 drjsanger@github
 ***********************************************************************/

#include <iostream>
#include <vector>
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

extern "C" {
  void add_new_tag(std::map<std::string, std::vector<int> > *dict, const char* new_key, int beg, int end)
  {
    // //cout<< "Added tag: " << new_key << endl;
    std::vector<int> counts((end-beg+1)*11*2,0);
    dict->insert({new_key, counts});
  }
  void bam2R_pileup_function(const bam_pileup1_t *pl, int pos, int n_plp, nttable_t& nttable,std::map<std::string,std::vector<int> >& dic)
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
	    int bc_absent;
	    const char *rg;
	    rg = (char *)bam_aux_get(p->b, "CB");	    
	    if ( !rg ) {
	      Rf_warning("Barcode is missing in some reads!");
	      missingBC++;
	    }else{
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
		if(p->qpos < nttable.head_clip || (bam_is_rev(p->b) && ((bam1_core_t)p->b->core).l_qseq - p->qpos < nttable.head_clip)){
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
  SEXP bam2R_10x(SEXP bamfiletest, SEXP reftest, SEXP begtest, SEXP endtest, int countstest[], SEXP qtest, SEXP mqtest, SEXP stest, SEXP head_cliptest, SEXP maxdepthtest, SEXP verbosetest, SEXP masktest, SEXP keepflagtest, SEXP maxmismatchestest)
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
    
    bam_plp_t buf = NULL;
    bam1_t *b = NULL;
    hts_itr_t *iter = NULL;
    bam_hdr_t *head = NULL;

    int c = 0;
    nttable_t nttable;
    nttable.q = *q; //Base quality cutoff
    nttable.s = *s; //Strand (2=both)
    nttable.head_clip = *head_clip;
    nttable.i = 0;
    
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
      Rf_warning("BAM/CRAM is missing NM tag, ignoring max.mismatches argument.\n");
      //return Rf_ScalarLogical(FALSE);
    }
    if (strcmp(ref, "") == 0) { // if a region is not specified
      //Replicate sampileup functionality (uses above mask without supplementary)
      while((ret = sam_read1(nttable.in, head, b)) >= 0){
  	if ((b->core.flag & *mask)==0 && b->core.qual >= *mq && (b->core.flag & *keepflag)==*keepflag && ((paux && bam_aux2i(bam_aux_get(b,"NM")) <= *maxmismatches ) || !paux || *maxmismatches == -1)){	 
	  bam_plp_push(buf, b);
  	};
  	while ( (pl=bam_plp_next(buf, &tid, &pos, &n_plp)) != 0) {
  	  bam2R_pileup_function(pl,pos,n_plp,nttable,dic);
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
	Rf_error("Invalid sequence %s\n", ref);
  	return Rf_ScalarLogical(FALSE);//return 1;
      }
      char *region = NULL;
      region = (char*) malloc(sizeof(ref)+sizeof(":")+sizeof("-")+(sizeof(char)*50));      
      sprintf(region,"%s:%d-%d",ref,nttable.beg,nttable.end);
      if(*verbose){
	Rprintf("Reading %s at coordinates %s.\n", bamfile, region);
      }
      //Implement a fetch style iterator
      hts_itr_t *iter = sam_itr_querys(idx, head, region);
      int result;
      while ((result = sam_itr_next(nttable.in, iter, b)) >= 0) {
  	if ((b->core.flag & *mask)==0 && b->core.qual >= *mq && (b->core.flag & *keepflag)==*keepflag && ((paux && bam_aux2i(bam_aux_get(b,"NM")) <= *maxmismatches ) || !paux || *maxmismatches == -1)){
  	  bam_plp_push(buf, b);
  	};
  	while ( (pl=bam_plp_next(buf, &tid, &pos, &n_plp)) != 0) {
  	  bam2R_pileup_function(pl,pos,n_plp,nttable,dic);
  	}
      }
      if(result < -1){
	Rf_error("Error code (%d) encountered reading sam iterator.\n", result);
  	return Rf_ScalarLogical(FALSE);//return 1;
      }
      free(region);
      sam_itr_destroy(iter);
      hts_idx_destroy(idx);
    }

    bam_plp_push(buf,0); // finalize pileup

    while ( (pl=bam_plp_next(buf, &tid, &pos, &n_plp)) != 0) {
      bam2R_pileup_function(pl,pos,n_plp,nttable, dic);
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
			     {"bam2R_10x", (DL_FUNC) &bam2R_10x, 12}
  };
  
  void R_init_bam2R_10x(DllInfo *info) {
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  }
} // extern "C"
