#ifndef _ND_H_
#define _ND_H_

/***************************************************
* File: nd.h
* This file consists of Nested Dissection 
*                       Method subroutines
***************************************************/

#include "common.h"
#include "core.h"

// define epsilon (small value)
#define ND_EPS_THRESHOLD           1

//  Nested Dissection
int  nd     (TMatrix_DCSR* /*matr*/, real /*threshold*/);
int  nd_mod (TMatrix_DCSR* /*matr*/);
int  nd_orig(TMatrix_DCSR* /*matr*/);

#endif
