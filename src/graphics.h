#ifndef _GRAPHICS_H_
#define _GRAPHICS_H_

/***************************************************
* File: graphics.h
* Set of function for making graphic portrain of
* matricies in PNG-format
***************************************************/

#include "core.h"

// define max size of the picure: N x N
#define MAX_PNG_SIZE 640

int make_matrix_portrait(TMatrix_DCSR* /*matr*/, const char* /*filename*/);

#endif
