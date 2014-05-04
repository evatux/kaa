#ifndef _CORE_H_
#define _CORE_H_

/***************************************************
* File: core.h
* Realization primitives:
*   o) matrix formats (CSR, DiagCSR)
*   o) graphs (w/wo weight)
***************************************************/

#include "common.h"

typedef struct {
    real *val;
    int  rows;
    int  cols;
} TMatrix_SMP;

typedef struct {
    real *val;
    int  *col_ind;
    int  *row_ptr;
    int  nonz;
    int  rows;
    int  cols;
} TMatrix_CSR;

typedef struct {
    real *val;
    int  size;
} TVector_SMP;

//  Matrix CSR
int  matrix_create      (TMatrix_CSR* /*matr*/, int /*rows*/, int /*cols*/, int /*nonz*/, int /*clean_flag*/);
int  matrix_load        (TMatrix_CSR* /*matr*/, const char* /*filename*/);
int  matrix_load_fmc    (TMatrix_CSR* /*matr*/, const char* /*filename*/);
int  matrix_save        (TMatrix_CSR* /*matr*/, const char* /*filename*/);
void matrix_destroy     (TMatrix_CSR* /*matr*/);

int  matrix_portrait    (TMatrix_CSR* /*matr*/, const char* /*filename*/);
void matrix_show        (TMatrix_CSR* /*matr*/, int /*flag_ordered*/);


//  Matrix SMP
int  matrix_smp_create  (TMatrix_SMP* /*matr*/, int /*rows*/, int /*cols*/);
void matrix_smp_destroy (TMatrix_SMP* /*matr*/);

void matrix_smp_show    (TMatrix_SMP* /*matr*/);


//  Vector_SMP
int  vector_create      (TVector_SMP* /*vect*/, int /*size*/);
void vector_show        (TVector_SMP* /*vect*/);
void vector_destroy     (TVector_SMP* /*vect*/);


//  Matrix-Matrix interface
int  matrix_copy        (TMatrix_CSR* /*src*/, TMatrix_CSR* /*dst*/);
int  matrix_perm_copy   (TMatrix_CSR* /*src*/, TMatrix_CSR* /*dst*/, int* /*perm_rows*/, int* /*perm_cols*/);
int  matrix_smp2csr     (TMatrix_SMP* /*src*/, TMatrix_CSR* /*dst*/);
int  matrix_csr2smp     (TMatrix_CSR* /*src*/, TMatrix_SMP* /*dst*/);


//  Matrix-Vector interface
int  matrix_vector_mul  (TMatrix_CSR* /*matr*/, TVector_SMP* /*in*/,TVector_SMP* /*out*/);

#endif
