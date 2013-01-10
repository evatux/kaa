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
int  find_periphery  (TWGraph* /*gr*/, int* /*per*/);
int  find_permutation(TWGraph* /*gr*/, int  /*root*/, int** /*_perm*/, int** /*_invp*/);
int  graph_reorder   (TWGraph* /*gr*/, int* /*perm*/, int* /*invp*/);

//	RCM
int  rcm(TMatrix_DCSR* /*matr*/);

#endif

