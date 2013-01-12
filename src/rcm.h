#ifndef _RCM_H_
#define _RCM_H_

/***************************************************
* File: rcm.h
* This file consists of RCM subroutines
***************************************************/

#include "common.h"
#include "core.h"

// define epsilon (small value)
#define EPS_THRESHOLD			1

//	RCM subroutines
int  find_permutation(TWGraph* /*gr*/, int** /*_perm*/, int** /*_invp*/, real /*threshold*/);
int  graph_reorder   (TWGraph* /*gr*/, int*  /* perm*/, int*  /* invp*/);

//	RCM
int  rcm     (TMatrix_DCSR* /*matr*/, real /*threshold*/);
int  rcm_mod (TMatrix_DCSR* /*matr*/);
int  rcm_orig(TMatrix_DCSR* /*matr*/);

#endif

