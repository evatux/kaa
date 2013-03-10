#ifndef _GRAPHICS_H_
#define _GRAPHICS_H_

/***************************************************
* File: graphics.h
* Set of function for making graphic portrain of
* matricies in PNG-format
***************************************************/

#include "core.h"

// define max size of the picure: N x N
#define MAX_PNG_SIZE 1024

#ifdef _TEST_DOUBLE
#define GRAPH_ZERO_THRESHOLD 1e-7
#else
#define GRAPH_ZERO_THRESHOLD 1e-4
#endif

// define colors (r,g,b)
#define CL_BLACK    0
#define CL_RED      4
#define CL_GREEN    2
#define CL_BLUE     1
#define CL_WHITE    7

int make_matrix_portrait(TMatrix_DCSR* /*matr*/, const char* /*filename*/);
int make_matrix_portrait_color(TMatrix_DCSR* /*matr*/, const char* /*filename*/, real /*threshold*/);

#endif
