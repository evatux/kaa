#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "core.h"
#include "spmv.h"
#include "bench.h"

static real foo(int i, int j, int rows, int cols, double sparsity)
{
    double x   = (sin(i) * cols / 1.387 + cos(j) * rows / 2.713);
    double phi = x - ((int)(x/M_PI/2))*2*M_PI;

    if (fabs(sin(phi)) > fabs(sin(sparsity*M_PI/2.))) return 0;

    return 3. + cos(i + cols) + sin(j + rows);
}

int correctness(spmv_kernel_t *ker, info_t info, corr_t *c)
{
    int rows, cols;

    TMatrix_CSR *matr;

    if (c->corr_type == CORRECTNESS_FAST) {
        rows = correctness_rows;
        cols = correctness_cols;
        double sparsity = correctness_sparsity;

        TMatrix_SMP matr_smp;
        DSAFE(matrix_smp_create(&matr_smp, rows, cols));

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                matr_smp.val[i*cols + j] =
                    foo(i, j, rows, cols, sparsity);
            }
        }

        TMatrix_CSR _matr;
        matr = &_matr;
        DSAFE(matrix_smp2csr(&matr_smp, matr));
        matrix_smp_destroy(&matr_smp);

        D(
            char debug_name[MAX_FILENAME_LENGTH];
            sprintf(debug_name, "debug_correctness_%d_%d_%d_%.2f.png",
                rows, cols, matr->nonz, sparsity);
            matrix_portrait(matr, debug_name);
        );
    } else
    {
        matr = c->matr;
        rows = matr->rows;
        cols = matr->cols;
    }

    TVector_SMP in, out, out_true;
    DSAFE(vector_create(&in,        cols));
    DSAFE(vector_create(&out,       rows));
    DSAFE(vector_create(&out_true,  rows));

    for (int i = 0; i < cols; ++i) in.val[i]  = foo(i, 0, cols, 0, 1.);
    for (int i = 0; i < rows; ++i) out.val[i] = out_true.val[i] = 0.;

    DSAFE(matrix_vector_mul(matr, &in, &out_true));

    desc_t *desc;
    DSAFE(ker->create(&desc, matr, info));
    DSAFE(ker->compute(desc, &in,  &out));
    ker->destroy(desc);
    if (c->corr_type == CORRECTNESS_FAST)
        matrix_destroy(matr);

    double err = 0.;
    for (int i = 0; i < rows; ++i) {
        double diff = abs(out.val[i] - out_true.val[i]);
        if (abs(out_true.val[i]) > EPS)
            diff /= abs(out_true.val[i]);
        if (diff > err) err = diff;
    }
    c->err  = err;
    c->fail = (err > THRESHOLD);
    return ERROR_NO_ERROR;
}
