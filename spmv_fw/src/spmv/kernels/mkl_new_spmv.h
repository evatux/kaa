#ifndef _MKL_SPMV_H_
#define _MKL_SPMV_H_

#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include "spmv.h"
#include "spmv_interface.h"

#ifdef _DOUBLE
#define mkl_cspblas_gemv    mkl_cspblas_dcsrgemv
#define sparseXcsr2csr      sparseDcsr2csr
#define sparseXcsrmv        sparseDcsrmv
#else
#define mkl_cspblas_gemv    mkl_cspblas_scsrgemv
#define sparseXcsr2csr      sparseScsr2csr
#define sparseXcsrmv        sparseScsrmv
#endif

typedef struct {
    int hash;
    sparseMatDescr_t    desc;
    sparseCSRMatrix_t   matr;
} mkl_new_desc_t;

#endif
