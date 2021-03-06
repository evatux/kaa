#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "common.h"
#include "core.h"
#include "spmv.h"
#include "ref.h"

#define REF_HASH 0x00ff00ff

static int  mv_omp_mul(TMatrix_CSR *matr, TVector_SMP *in, TVector_SMP *out)
{
    int rows     = matr->rows;
    int cols     = matr->cols;
    int in_size  = in->size;
    int out_size = out->size;

    #pragma omp parallel for
    for (int i = 0; i < rows; ++i) {
        real s = 0.;
        for (int j = matr->row_ptr[i]; j < matr->row_ptr[i+1]; ++j)
        {
            s += matr->val[j] * in->val[matr->col_ind[j]];
        }
        out->val[i] = s;
    }

    return ERROR_NO_ERROR;
}

static int ref_desc_create(ref_desc_t **_desc, TMatrix_CSR *matr, info_t info)
{
    ref_desc_t *desc = malloc(sizeof(ref_desc_t));
    if (desc == NULL) return DE(ERROR_MEMORY);
    int err = matrix_copy(matr, &desc->matr);
    if (err) {
        *_desc = NULL;
        free(desc);
        return DE(err);
    }
    desc->hash = REF_HASH;
    *_desc = desc;
    return ERROR_NO_ERROR;
}

static int ref_compute(ref_desc_t *desc, TVector_SMP *in, TVector_SMP *out)
{
    if (desc->hash != REF_HASH) return DE(ERROR_INVALID_DESC);

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

    DSAFE(mv_omp_mul(matr, in, out));

    return ERROR_NO_ERROR;
}

static void ref_desc_destroy(ref_desc_t *desc)
{
    if (desc == NULL) return;
    if (desc->hash != REF_HASH) return;

    matrix_destroy(&desc->matr);
    free(desc);
}

spmv_kernel_t ref_omp()
{
    spmv_kernel_t ker = {
        (spmv_desc_create_t)ref_desc_create,
        (spmv_compute_t)ref_compute,
        (spmv_desc_destroy_t)ref_desc_destroy,
        "ref_omp_code"
    };
    return ker;
}
