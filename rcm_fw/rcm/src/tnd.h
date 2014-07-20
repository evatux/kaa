#ifndef _TND_H_
#define _TND_H_

/***************************************************
* File: nd.h
* This file consists of Nested Dissection 
*                       Method subroutines
***************************************************/

#include "common.h"
#include "core.h"

// define epsilon (small value)
#define TND_EPS_THRESHOLD           0

#if (TND_SMART_TYPE != 1 && TND_SMART_TYPE != 2)
#define TND_SMART_TYPE 2
#endif

//  Nested Dissection
int  tnd     (TMatrix_DCSR* /*matr*/, real /*threshold*/);
int  tnd_perm(TWGraph* /*gr*/, int** /*_perm*/, int** /*_invp*/);
int  tnd_mod (TMatrix_DCSR* /*matr*/);
int  tnd_orig(TMatrix_DCSR* /*matr*/);

#endif
