#ifndef _RMD_H_
#define _RMD_H_

/***************************************************
* File: rmd.h
* This file consists of Recursive Minimal Degree subroutines
***************************************************/

#include "common.h"
#include "core.h"

// define epsilon (small value)
#define RMD_EPS_THRESHOLD           1

//  RMD
int rmd     (TMatrix_DCSR* /*matr*/, real /*threshold*/);
int rmd_mod (TMatrix_DCSR* /*matr*/);
int rmd_orig(TMatrix_DCSR* /*matr*/);

#endif
