#include "solver.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "common.h"
#include "core.h"

static real vector_normalize(real *X, int size)
{
    int i;
    real norm = 0;

#ifdef _OPENMP
#pragma omp parallel for reduction (+:norm)
#endif
    for(i = 0; i < size; ++i) {
        norm += X[i]*X[i];
    }

    norm = SQRT(norm);
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
    for(i = 0; i < size; ++i) {
        X[i] /= norm;
    }

    return norm;
}

// Y <- AX
int matrix_vector_mult(TMatrix_DCSR *A, real *X, real *Y)
{
    int i, j;
#ifdef _OPENMP
#pragma omp parallel for private(i,j)
#endif
    for (i = 0; i < A->size; ++i) {
        Y[i] = A->diag[i]*X[i];
        for (j = A->row_ptr[i]; j < A->row_ptr[i+1]; j++)
            Y[i] += A->val[j] * X[A->col_ind[j]];
    }
    return ERROR_NO_ERROR;
}

// Solve: Y <- (A^-1)Y, cholesky(A) = LD
int solver(TMatrix_DCSR *LD, real *Y)
{
    int i, j, dbz = 0;

    // Y <- (L^-1)Y
    for (i = 0; i < LD->size; ++i) {
        for (j = LD->row_ptr[i]; LD->col_ind[j] < i && j < LD->row_ptr[i+1]; ++j)
            Y[i] -= LD->val[j] * Y[LD->col_ind[j]];
    }

    // Y <- (D^-1)Y
#ifdef _OPENMP
#pragma omp parallel for private(i) shared(dbz)
#endif
    for (i = 0; i < LD->size; ++i) 
        if ( FABS(LD->diag[i]) < INV_EPS ) dbz = 1;
        else Y[i] /= LD->diag[i];
    if (dbz == 1) {
#ifdef _DEBUG_LEVEL_SOLVER
        printf("[debug] {solver}: division by zero found\n");
#endif
        return ERROR_DIV_BY_ZERO;
    }

    // Y <- (L^-T)Y
    for (i = LD->size-1; i >= 0; --i) {
        for (j = LD->row_ptr[i]; LD->col_ind[j] < i && j < LD->row_ptr[i+1]; ++j) 
            Y[LD->col_ind[j]] -= LD->val[j] * Y[i];
    }

    return ERROR_NO_ERROR;
}

// Make: E` <- A^-1 A, where A^-1 is made from cholesky of A
int make_ident(TMatrix_DCSR *A, TMatrix_DCSR *LD, TMatrix_DCSR *E)
{
#define CLEANUP() do { if (y) free(y); if (ES.val) free(ES.val); } while(0)
    int size = A->size;
    int i, j, it;
    int err;
    real *y;
    TMatrix_Simple ES;

    y = (real*)malloc(sizeof(real)*size);
    ES.val = (real*)malloc(sizeof(real)*size*size);
    ES.size = size;

//matrix_show(LD, 1);

    if (y == NULL || ES.val == NULL) {
        CLEANUP();
        fprintf(stderr, "error [make_ident]: memory allocation error\n");
        return ERROR_MEMORY_ALLOCATION;
    }

    // y <- A[*, j], y <- (A^-1) y
    for (j = 0; j < size; ++j)
    {
        for (i = 0; i < size; ++i) y[i] = 0;
        y[j] = A->diag[j];
        for (i = 0; i < size; ++i)
        {
            for (it = A->row_ptr[i]; it < A->row_ptr[i+1]; ++it)
            {
                if (A->col_ind[it] == j) {
                    y[i] = A->val[it];
//if (FABS(y[i]) > 1e-3) printf("%.2e ", y[i]);
                    break;
                }
            }
        }

        solver(LD, y);

        for (i = 0; i < size; ++i) ES.val[i*size + j] = y[i];
    }

    err = matrix_convert_simp2dcsr(&ES, E);
    if (err != ERROR_NO_ERROR) {
        fprintf(stderr, "error [make_ident]: matrix convert error (%d)\n", err);
        PRINT_ERROR_MESSAGE(err);
    }

    CLEANUP();
    return err;
}

