#ifndef _SPMV_H_
#define _SPMV_H_

#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "core.h"

/*
 * SpMV API
 */

int  (*spmv_desc_create_t) (spmv_desc_t *, TMatrix_CSR *, void *);
int  (*spmv_compute_t)     (spmv_desc_t *, TVector_SMP *, TVector_SMP *);
void (*spmv_desc_destroy_t)(spmv_desc_t *);

typedef struct {
    spmv_desc_create_t  create;
    spmv_compute_t      compute;
    spmv_desc_destroy_t destroy;
} spmv_kernel_t;

#endif
