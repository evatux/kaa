#ifndef _SPMV_H_
#define _SPMV_H_

#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "core.h"

#define MAX_INFO_SIZE   128
#define MAX_KER_NUM     10
#define KER_NONE        (-1)

//  SpMV API

typedef void *desc_t;
typedef char info_t[MAX_INFO_SIZE];

typedef int  (*spmv_desc_create_t) (desc_t **, TMatrix_CSR *, info_t);
typedef int  (*spmv_compute_t)     (desc_t *,  TVector_SMP *, TVector_SMP *);
typedef void (*spmv_desc_destroy_t)(desc_t *);

typedef struct {
    spmv_desc_create_t  create;
    spmv_compute_t      compute;
    spmv_desc_destroy_t destroy;
    const char*         name;
} spmv_kernel_t;

typedef spmv_kernel_t (*spmv_kernel_foo_t)();


//  SpMV kernels

spmv_kernel_t ref();
spmv_kernel_t mkl_spmv();

static spmv_kernel_foo_t kernel_list[] = {
    ref,
#ifdef _MKL_ENABLE
    mkl_spmv,
#endif
};

extern const int kernel_num;


//  SpMV utils

double gflop(TMatrix_CSR* /*matr*/);

#endif