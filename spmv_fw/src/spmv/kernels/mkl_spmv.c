#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include "common.h"
#include "core.h"
#include "spmv.h"
#include "mkl_spmv.h"

#define MKL_SPMV_HASH 0x11515111

static int mkl_desc_create(mkl_desc_t **_desc, TMatrix_CSR *matr, info_t info)
{
    mkl_desc_t *desc = malloc(sizeof(mkl_desc_t));
    if (desc == NULL) return DE(ERROR_MEMORY);
    int err = matrix_copy(matr, &desc->matr);
    if (err) {
        *_desc = NULL;
        free(desc);
        return DE(err);
    }
    desc->hash = MKL_SPMV_HASH;
    *_desc = desc;
    return ERROR_NO_ERROR;
}

static int mkl_compute(mkl_desc_t *desc, TVector_SMP *in, TVector_SMP *out)
{
    if (desc->hash != MKL_SPMV_HASH) return DE(ERROR_INVALID_DESC);

    TMatrix_CSR *matr = &desc->matr;
    int rows     = matr->rows;
    int cols     = matr->cols;
    int in_size  = in->size;
    int out_size = out->size;

    if (!(rows == out_size && cols == in_size))
    {
        D(fprintf(stderr, "matr: %d x %d\nin: %d, out: %d\n",
                    rows, cols, in_size, out_size));
        return DE(ERROR_INVALID_CONF);
    }

    MKL_INT m = rows;
    MKL_INT *ia = (MKL_INT*)matr->row_ptr;
    MKL_INT *ja = (MKL_INT*)matr->col_ind;
    mkl_cspblas_gemv("N", &m, matr->val, ia, ja, in->val, out->val);

    return ERROR_NO_ERROR;
}

static void mkl_desc_destroy(mkl_desc_t *desc)
{
    if (desc == NULL) return;
    if (desc->hash != MKL_SPMV_HASH) return;

    matrix_destroy(&desc->matr);
    free(desc);
}

spmv_kernel_t mkl_spmv()
{
    spmv_kernel_t ker = {
        (spmv_desc_create_t)mkl_desc_create,
        (spmv_compute_t)mkl_compute,
        (spmv_desc_destroy_t)mkl_desc_destroy,
        "mkl_spmv"
    };
    return ker;
}
