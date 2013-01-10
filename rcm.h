#ifndef _RCM_H_
#define _RCM_H_

#ifdef  _DEBUG
#define _DEBUG_LEVEL_0
#define _DEBUG_LEVEL_1
#endif

#define PRINT_ZEROES

/*
 * Notation: 
 *		Numeration starts from 0 to n-1
 *
 */

//#define _TEST_DOUBLE 

#ifdef _TEST_DOUBLE
typedef double real;
#define FABS fabs
#else
typedef float real;
#define FABS fabsf
#endif

#define ERROR_NO_ERROR			0
#define ERROR_MEMORY_ALLOCATION 	2
#define ERROR_FILE_IO 			17

#define EPS_THRESHOLD			1

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

//	RCM subroutines
int  find_periphery(TWGraph* /*gr*/, int* /*per*/);
int  find_permutation(TWGraph* /*gr*/, int /*root*/, int** /*_perm*/, int** /*_invp*/);
int  graph_reorder(TWGraph* /*gr*/, int* /*perm*/, int* /*invp*/);

//	RCM
int  rcm(TMatrix_DCRS* /*matr*/);

//	old interface
void old_interface_matrix_load(TMatrix_CRS* /*matr*/, const char* /*filename*/);
void old_interface_matrix_show(TMatrix_CRS* /*matr*/);
void old_interface_graph_builder(TGraph* /*gr*/, TMatrix_CRS* /*matr*/);
int  old_interface_find_periphery(TGraph* /*gr*/);
void old_interface_find_permutation(TGraph* /*gr*/, int /*root*/, int* /*perm*/);
void old_interface_rcm(TMatrix_CRS* /*matr*/);

#endif

