#ifndef _CORE_H_
#define _CORE_H_

/***************************************************
* File: core.h
* Realization primitives:
*	o) matrix formats (CSR, DiagCSR)
*	o) graphs (w/wo weight)
***************************************************/

#include "common.h"

typedef struct {
	int *adjncy;
	int *xadj;
	int size;
} TGraph;

typedef struct {
	real *wedge;
	int  *adjncy;
	int  *xadj;
	real *wvert;
	int  size;
	int  nonz;
} TWGraph;

typedef struct {
	real *val;
	int  *col_ind;
	int  *row_ptr;
	int  nonz;
	int  size;
} TMatrix_CSR;

typedef struct {
	real *diag;

	real *val;
	int  *col_ind;
	int  *row_ptr;
	
	int  nonz;
	int  size;
} TMatrix_DCSR;

//	Matrix & Graph interface
int  matrix_load(TMatrix_DCSR* /*matr*/, const char* /*filename*/);
int  matrix_save(TMatrix_DCSR* /*matr*/, const char* /*filename*/);
int  matrix_portrait(TMatrix_DCSR* /*matr*/, const char* /*filename*/, real /*threshold*/);
int  matrix_portrait_pattern(TMatrix_DCSR* /*matr*/, const char* /*pattern*/, const char* /*suffix*/, real /*threshold*/);

void matrix_show(TMatrix_DCSR* /*matr*/, int /*flag_ordered*/);
void graph_show (TWGraph* /*gr*/);

//	Matrix <--> Graph interface
int  build_graph (TWGraph* /*gr*/, TMatrix_DCSR* /*matr*/);
int  build_matrix(TWGraph* /*gr*/, TMatrix_DCSR* /*matr*/, int /*flag_new*/);

#endif
