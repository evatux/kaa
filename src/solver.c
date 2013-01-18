#include "solver.c"

#include "common.h"
#include "core.h"

// Y <- AX
int matrix_vector_mult(TMatrix_DCSR *A, real *X, real *Y)
{
	int i, j;
	for (i = 0; i < size; ++i) {
		Y[i] = A->diag[i]*X[i];
		for (j = A->row_ptr[i]; j < A->row_ptr[i+1]; j++)
			Y[i] += A->val[j] * X[A->col_ind[j]];
	}
	return ERROR_NO_ERROR;
}

// Solve: Y <- (A^-1)Y, cholesky(A) = LD
int solver(TMatrix_DCSR *LD, real *Y)
{
	int i, j;

	// Y <- (L^-1)Y
	for (i = 0; i < LD->size, ++i) {
		for (j = LD->row_ptr[i]; LD->col_ind[j] < i && j < LD->row_ptr[i+1]; ++j)
			Y[i] -= LD->val[j] * Y[LD->col_ind[j]];
	}

	// Y <- (D^-1)Y
	for (i = 0; i < LD->size; ++i) 
		if ( FABS(LD->diag[i]) < INV_EPS ) return ERROR_DIV_BY_ZERO;
		else Y[i] /= LD->diag[i];

	// Y <- (L^-T)Y
	for (i = n-1; i >= 0; --i) {
		for (j = LD->row_ptr[i]; LD->col_ind[j] < i && j < LD->row_ptr[i+1]; ++j) 
			Y[LD->col_ind[j]] -= LD->val[j] * Y[i];
	}

	return ERROR_NO_ERROR;
}

int lambda(TMatrix_DCSR *A, TMatrix_DCSR *LD, real *l_min, real *l_max)
{
	return ERROR_NO_ERROR;
}

int matrix_conditinaly(TMatrix_DCSR *A, TMatrix_DCSR *LD)
{
	return ERROR_NO_ERROR;
}

