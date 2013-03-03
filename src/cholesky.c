#include "cholesky.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

#include "core.h"
#include "common.h"

#define EXCEPTION(cond, msg, ind, EXCP) \
    do {                                \
        if ((cond)) {                   \
            fprintf(stderr, msg, ind);  \
            PRINT_ERROR_MESSAGE(EXCP);  \
            exit(EXCP);                 \
        }                               \
    } while(0)

int cholesky_decomposition(TMatrix_DCSR *A, TMatrix_DCSR *LD, const real cheps, int *neps)
{
    int err;
    int i, j, k;
    int size = A->size;
    real sum, cur, ljj;
    real* LS;
    TMatrix_Simple L_simp;

    *neps = 0;  // number of convertations 0 --> cheps

//matrix_show(A, 0);
    err = matrix_convert_dcsr2simp(A, &L_simp);
    if (err != ERROR_NO_ERROR) return err;
//matrix_simp_show(&L_simp);
    LS = L_simp.val;

    // Cholesky decomposition
    for (j = 0; j < size; ++j) {
        ljj = LS[j*size+j]; for (k=0; k<j; ++k) ljj -= LS[j*size+k]*LS[j*size+k]*LS[k*size+k];
//EXCEPTION(sum < 0, "l[j,j] < 0 ; j = %d\n", j, ERROR_NEGATIVE_SQRT);

        if (FABS(ljj) < cheps) {    // this is modification of cholesky
            ljj = LS[j*size+j] = cheps*SGN(ljj);
            *neps += 1;
        } else LS[j*size+j] = ljj;

        for (i = j+1; i < size; ++i) {
            sum = 0.; for (k=0; k<j; ++k) sum += LS[i*size+k]*LS[j*size+k]*LS[k*size+k];
            LS[i*size+j] = (LS[i*size+j] - sum) / ljj; 
        }
    }

//matrix_simp_show(&L_simp);
    for (i = 0; i < size; ++i)
        for (j = i + 1; j < size; ++j)
            LS[i*size+j] = 0.;

    err = matrix_convert_simp2dcsr(&L_simp, LD);
    if (err != ERROR_NO_ERROR) {
        matrix_simp_destroy(&L_simp);
        return err;
    }

    matrix_simp_destroy(&L_simp);
    return ERROR_NO_ERROR;
}
