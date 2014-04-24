#ifndef _WND_H_
#define _WND_H_

/***************************************************
* File: wnd.h
* This file consists of Weight Nested Dissection 
*                       Method subroutines
***************************************************/

#include "common.h"
#include "core.h"

// define epsilon (small value)
#define WND_EPS_THRESHOLD    1
#define WND_AM_THRESHOLD    10

//  Weight Nested Dissection
int  wnd     (TMatrix_DCSR* /*matr*/, real /*threshold*/);
int  wnd_mod (TMatrix_DCSR* /*matr*/);
int  wnd_orig(TMatrix_DCSR* /*matr*/);

#endif
