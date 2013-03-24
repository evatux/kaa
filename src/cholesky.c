#include "cholesky.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

#include "core.h"
#include "unit.h"
#include "common.h"

#define EXCEPTION(cond, msg, ind, EXCP) \
    do {                                \
        if ((cond)) {                   \
            fprintf(stderr, msg, ind);  \
            PRINT_ERROR_MESSAGE(EXCP);  \
            exit(EXCP);                 \
        }                               \
    } while(0)

void cholesky_preanalysis(TMatrix_DCSR *A, real threshold, int *pre_neps)
{
    int res = 0;
    int ci, i, flag;
    for (i = 0; i < A->size; ++i) {
        if (FABS(A->diag[i]) < threshold) {
            flag = 1;
            for (ci = A->row_ptr[i]; ci < A->row_ptr[i+1]; ci++)
                if ((A->col_ind[ci] < i) && (A->val[ci] != 0.))
                {
                    flag = 0; 
                    break;
                }
            if (flag)
            {
                res++;
            }
        }
    }
    *pre_neps = res;
}

int cholesky_decomposition(TMatrix_DCSR *A, TMatrix_DCSR *LD, const real threshold, const real substitute, int *neps, int **neps_list)
{
    int err;
    int i, j, k;
    int size = A->size;
    real sum, cur, ljj;
    real* LS;
    TQueue queue;
    TMatrix_Simple L_simp;

    *neps = 0;  // number of convertations 0 --> cheps
    queue_init(&queue);

//matrix_show(A, 0);
    err = matrix_convert_dcsr2simp(A, &L_simp);
    if (err != ERROR_NO_ERROR) return err;
//matrix_simp_show(&L_simp);
    LS = L_simp.val;

    // Cholesky decomposition
    for (j = 0; j < size; ++j) {
#ifdef _DEBUG_LEVEL_CHOLESKY
printf("%.2e: ", LS[j*size+j]);
#endif
        ljj = LS[j*size+j]; for (k=0; k<j; ++k) {
#ifdef _DEBUG_LEVEL_CHOLESKY
            if (LS[j*size+k]*LS[j*size+k]*LS[k*size+k] != 0.) printf("%.2e ", LS[j*size+k]*LS[j*size+k]*LS[k*size+k]);
#endif
            ljj -= LS[j*size+k]*LS[j*size+k]*LS[k*size+k];
        }
//EXCEPTION(sum < 0, "l[j,j] < 0 ; j = %d\n", j, ERROR_NEGATIVE_SQRT);

        if (FABS(ljj) < threshold) {    // this is modification of cholesky
            ljj = LS[j*size+j] = substitute*SGN(ljj);
            *neps += 1;
            if (neps_list) queue_push(&queue, j);
        } else LS[j*size+j] = ljj;
#ifdef _DEBUG_LEVEL_CHOLESKY
printf(" --> %.2e\n", LS[j*size+j]);
#endif
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

    if ((*neps) != 0) {
        *neps_list = (int*)malloc(sizeof(int)*(*neps));
        if ((*neps_list) == NULL) {
            fprintf(stderr, "[cholesky] memory allocation error\n");
            return ERROR_MEMORY_ALLOCATION;
        }
        for (i = 0; i < *neps; ++i) {
            queue_pop(&queue, &j);
            (*neps_list)[i] = j;
        }
    } else *neps_list = NULL;

    return ERROR_NO_ERROR;
}
