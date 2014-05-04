#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <omp.h>
#include "common.h"
#include "core.h"
#include "spmv.h"
#include "gblock_omp.h"

#define GBLOCK_OMP_HASH 0x004B10C3

#define UNSET   (-2)

static double gblock_alpha = 0.1;

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

static int get_row_perm(TMatrix_CSR *matr, int **_row_perm)
{
    int err;
    int rows = matr->rows;
    int cols = matr->cols;
    int nonz = matr->nonz;

    real alpha = gblock_alpha;
    real cur_alpha;

    int  *perm      = (int* )malloc(sizeof(int )*rows);
    real *w         = (real*)malloc(sizeof(real)*cols);
    if (!perm || !wc_ind || !w) err = ERROR_MEMORY, goto fail;

    //  init perm
    for (int k = 0; k < rows; ++k) perm[k] = UNSET;

    //  init weight
    for (int j = 0; j < cols; ++j) w[j] = 0.;
    for (int c = 0; c < matr->row_ptr[1]; ++c)
        w[matr->col_ind[c]] = 1./matr->row_ptr[1];

    //  compute perm
    for (int k = 1; k < rows; ++k)
    {
        real bcorr = 0;
        int  brow = 0;

        for (int i = 1; i < rows; ++i)
        {
            //  skip, if already was here
            if (perm[i] != UNSET) continue;

            real corr = 0.;
            for (int c = matr->row_ptr[i]; c < matr->row_ptr[i+1]; ++c)
                corr += w[matr->col_ind[c]];
            if (corr > bcorr) {
                bcorr = corr;
                brow  = i;
            }
        }

        perm[brow] = k;

        if (corr > corr_threshold) {
            //  keeping history
            cur_alpha = alpha;
        } else {
            //  time to forget
            cur_alpha = 0.;
        }

        int l = matr->row_ptr[brow+1] - matr->row_ptr[brow];
        for (int j = 0; j < cols; ++j)
            w[j] *= cur_alpha;
        for (int c = matr->row_ptr[brow]; c < matr->row_ptr[brow+1]; ++c)
            w[matr->col_ind[c] - matr->row_ptr[blow]] += (1-cur_alpha)/l;
    }

    *_row_perm = perm;
    return ERROR_NO_ERROR;

fail:
    if (perm)   free(row_perm);
    if (wr)     free(wr);
    if (w)      free(w);
    *_row_perm = NULL;

    return DE(err);
}

static int gblock_desc_create(gblock_desc_t **_desc, TMatrix_CSR *matr, info_t info)
{
    int err;
    int *row_perm = NULL;
    int *col_perm = NULL;

    gblock_desc_t *desc = malloc(sizeof(gblock_desc_t));
    if (desc == NULL) return DE(ERROR_MEMORY);

    err = get_row_perm(matr, &row_perm);
    if (err) goto fail;

    err = get_col_perm(matr, &col_perm);
    if (err) goto fail;

    err = matrix_copy_with_perm(matr, &desc->matr, row_perm, col_perm);
    if (err) goto fail;

    desc->hash = GBLOCK_OMP_HASH;
    *_desc = desc;

    return ERROR_NO_ERROR;

fail:
    if (row_perm)   free(row_perm);
    if (col_perm)   free(col_perm);
    if (desc)       free(desc);
    *_desc = NULL;

    return DE(err);
}

static int gblock_compute(gblock_desc_t *desc, TVector_SMP *in, TVector_SMP *out)
{
    if (desc->hash != GBLOCK_OMP_HASH) return DE(ERROR_INVALID_DESC);

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

static void gblock_desc_destroy(gblock_desc_t *desc)
{
    if (desc == NULL) return;
    if (desc->hash != GBLOCK_OMP_HASH) return;

    matrix_destroy(&desc->matr);
    free(desc);
}

spmv_kernel_t gblock_omp()
{
    spmv_kernel_t ker = {
        (spmv_desc_create_t)gblock_desc_create,
        (spmv_compute_t)gblock_compute,
        (spmv_desc_destroy_t)gblock_desc_destroy,
        "gblock_omp_code"
    };
    return ker;
}
