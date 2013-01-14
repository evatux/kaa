#include "cholesky.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

#include "core.h"
#include "common.h"

#define EXCEPTION(cond, msg, ind, EXCP) \
	do { 								\
		if ((cond)) {					\
			fprintf(stderr, msg, ind);	\
			PRINT_ERROR_MESSAGE(EXCP);	\
			exit(EXCP);					\
		} 								\
	} while(0)

int cholesky_decomposition(TMatrix_DCSR *A, TMatrix_DCSR *L, const real cheps, int *neps)
{
	int err;
	int i, j, k;
	int size = A->size;
	real sum, cur, ljj;
	real* LS;
	TMatrix_Simple L_simp;

	*neps = 0;	// number of convertations 0 --> cheps

//matrix_show(A, 0);
	err = matrix_convert_dcsr2simp(A, &L_simp);
	if (err != ERROR_NO_ERROR) return err;
//matrix_simp_show(&L_simp);
	LS = L_simp.val;

	// Cholesky decomposition
	for (j = 0; j < size; ++j) {
		sum = LS[j*size+j]; for (k=0; k<j-1; ++k) sum -= LS[j*size+k]*LS[j*size+k];
EXCEPTION(sum < 0, "l[j,j] < 0 ; j = %d\n", j, ERROR_NEGATIVE_SQRT);

		ljj = LS[j*size+j] = SQRT(sum);

		if (ljj < cheps) {	// this is modification of cholesky
			ljj = LS[j*size+j] = cheps;
			*neps += 1;
		}

		for (i = j+1; i < size; ++i) {
			sum = 0.; for (k=0; k<j-1; ++k) sum += LS[i*size+k]*LS[j*size+k];
			LS[i*size+j] = (LS[i*size+j] - sum) / ljj; 
		}
	}

//matrix_simp_show(&L_simp);
	for (i = 0; i < size; ++i)
		for (j = i + 1; j < size; ++j)
			LS[i*size+j] = 0.;

	err = matrix_convert_simp2dcsr(&L_simp, L);
	if (err != ERROR_NO_ERROR) {
		matrix_simp_destroy(&L_simp);
		return err;
	}

	matrix_simp_destroy(&L_simp);
	return ERROR_NO_ERROR;
}
