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
} TMatrix_CRS;

typedef struct {
	real *diag;

	real *val;
	int  *col_ind;
	int  *row_ptr;
	
	int  nonz;
	int  size;
} TMatrix_DCRS;

//	Matrix & Graph interface
int  matrix_load(TMatrix_DCRS* /*matr*/, const char* /*filename*/);
int  matrix_save(TMatrix_DCRS* /*matr*/, const char* /*filename*/);
void matrix_show(TMatrix_DCRS* /*matr*/, int /*flag_ordered*/);
void graph_show (TWGraph* /*gr*/);

//	Matrix <--> Graph interface
int  build_graph (TWGraph* /*gr*/, TMatrix_DCRS* /*matr*/);
int  build_matrix(TWGraph* /*gr*/, TMatrix_DCRS* /*matr*/, int /*flag_new*/);

#endif
