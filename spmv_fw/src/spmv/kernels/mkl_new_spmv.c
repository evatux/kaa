#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include "common.h"
#include "core.h"
#include "spmv.h"
#include "spmv_interface.h"
#include "mkl_new_spmv.h"

#define MKL_NEW_SPMV_HASH 0x11514E11

static int mkl_desc_create(mkl_new_desc_t **_desc, TMatrix_CSR *matr, info_t info)
{
    mkl_new_desc_t *desc = malloc(sizeof(mkl_new_desc_t));
    if (desc == NULL) return DE(ERROR_MEMORY);

    sparseStatus_t status;
    status = sparseCreateCSRMatrix(&desc->matr, SPARSE_SCHEDULE_STATIC);
    if (status) goto fail;
    status = sparseCreateMatDescr(&desc->desc);
    if (status) goto fail;
    status = sparseXcsr2csr ( matr->rows, matr->cols, desc->desc, matr->val, matr->row_ptr, matr->col_ind, desc->matr );
    if (status) goto fail;

    desc->hash = MKL_NEW_SPMV_HASH;
    *_desc = desc;
    return ERROR_NO_ERROR;

fail:
    *_desc = NULL;
    free(desc);
    return DE(status);
}

static int mkl_compute(mkl_new_desc_t *desc, TVector_SMP *in, TVector_SMP *out)
{
    if (desc->hash != MKL_NEW_SPMV_HASH) return DE(ERROR_INVALID_DESC);

    int in_size  = in->size;
    int out_size = out->size;

    /*
    if (!(rows == out_size && cols == in_size))
    {
        D(fprintf(stderr, "matr: %d x %d\nin: %d, out: %d\n",
                    rows, cols, in_size, out_size));
        return DE(ERROR_INVALID_CONF);
    }
    */

    real alpha = 1;
    real beta  = 0;
    sparseStatus_t status;
    status = sparseXcsrmv (SPARSE_OPERATION_NON_TRANSPOSE,
            &alpha, desc->matr, in->val, &beta, out->val);

    if (status) return DE(1);
    return ERROR_NO_ERROR;
}

static void mkl_desc_destroy(mkl_new_desc_t *desc)
{
    if (desc == NULL) return;
    if (desc->hash != MKL_NEW_SPMV_HASH) return;

    sparseDestroyCSRMatrix (desc->matr);
    free(desc);
}

spmv_kernel_t mkl_new_spmv()
{
    spmv_kernel_t ker = {
        (spmv_desc_create_t)mkl_desc_create,
        (spmv_compute_t)mkl_compute,
        (spmv_desc_destroy_t)mkl_desc_destroy,
        "mkl_new_spmv"
    };
    return ker;
}
