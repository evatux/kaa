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

//  Matrix & Graph interface
int  matrix_create      (TMatrix_CSR* /*matr*/, int /*rows*/, int /*cols*/, int /*nonz*/, int /*clean_flag*/);
int  matrix_load        (TMatrix_CSR* /*matr*/, const char* /*filename*/);
int  matrix_load_fmc    (TMatrix_CSR* /*matr*/, const char* /*filename*/);
int  matrix_save        (TMatrix_CSR* /*matr*/, const char* /*filename*/);

int  matrix_portrait    (TMatrix_CSR* /*matr*/, const char* /*filename*/);

void matrix_smp_show    (TMatrix_SMP* /*matr*/);
void matrix_show        (TMatrix_CSR* /*matr*/, int /*flag_ordered*/);

void matrix_smp_destroy (TMatrix_SMP* /*matr*/);
void matrix_destroy     (TMatrix_CSR* /*matr*/);

//  Matrix <--> Matrix interface
int  matrix_copy        (TMatrix_CSR* /*src*/, TMatrix_CSR* /*dst*/);
int  matrix_convert_smp2csr(TMatrix_SMP* /*src*/, TMatrix_CSR* /*dst*/);
int  matrix_convert_csr2smp(TMatrix_CSR* /*src*/, TMatrix_SMP* /*dst*/);

#endif
