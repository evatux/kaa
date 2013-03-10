#ifndef _MD_H_
#define _MD_H_

/***************************************************
* File: md.h
* This file consists of Minimal Degree subroutines
***************************************************/

#include "common.h"
#include "core.h"

// define epsilon (small value)
#define MD_EPS_THRESHOLD           1

//  MD
int  md     (TMatrix_DCSR* /*matr*/, real /*threshold*/);
int  md_mod (TMatrix_DCSR* /*matr*/);
int  md_orig(TMatrix_DCSR* /*matr*/);

#endif
