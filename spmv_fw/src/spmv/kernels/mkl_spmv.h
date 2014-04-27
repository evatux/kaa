#ifndef _MKL_SPMV_H_
#define _MKL_SPMV_H_

#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include "spmv.h"

#ifdef _DOUBLE
#define mkl_cspblas_gemv    mkl_cspblas_dcsrgemv
#else
#define mkl_cspblas_gemv    mkl_cspblas_scsrgemv
#endif

typedef struct {
    int hash;
    TMatrix_CSR matr;
} mkl_desc_t;

#endif
