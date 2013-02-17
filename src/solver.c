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

int lambda(TMatrix_DCSR *A, TMatrix_DCSR *LD, real *l_min, real *l_max)
{
#define CLEANUP() do { if (x) free(x); if (y) free(y); } while(0)
	int size = A->size;
	int i, iter;
	real *x, *y;
	real norm0, norm1, delta;
	real M, m;

	x = (real*)malloc(sizeof(real)*size);
	y = (real*)malloc(sizeof(real)*size);

	if (x == NULL || y == NULL) {
		CLEANUP();
		fprintf(stderr, "error [lambda]: memory allocation error\n");
		return ERROR_MEMORY_ALLOCATION;
	}

#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
	for (i = 0; i < size; ++i) x[i] = 3*sin(0.3*(i+1));
	norm1 = vector_normalize(x, size);

	iter = 0;
	do {
		iter++;

		norm0 = norm1;
		matrix_vector_mult(A, x, y);
		solver(LD, y);
#ifdef _DEBUG_LEVEL_SOLEVER_2
printf("[[debug]] (iter: %d) x[0..2] = %.2e; %.2e; %.2e\n", iter, x[0], x[1], x[2]);
printf("[[debug]] (iter: %d) y[0..2] = %.2e; %.2e; %.2e\n", iter, y[0], y[1], y[2]);
#endif
		SWAP(real*, x, y);
		norm1 = vector_normalize(x, size);
#ifdef _DEBUG_LEVEL_SOLVER_2
printf("[[debug]] (iter: %d) norm0 = %.2e, norm1 = %.2e\n", iter, norm0, norm1);
#endif

		delta = FABS(norm1 - norm0);
	} while( (iter < MAX_ITER) && (delta > LAMBDA_MAX_DELTA) );
	*l_max = M = norm1;
#ifdef _DEBUG_LEVEL_SOLVER
	printf("[debug] {solver}: l_max = %.2e [iter: %d]\n", M, iter);
#endif

#undef CLEANUP
	return ERROR_NO_ERROR;
}

int matrix_conditinaly(TMatrix_DCSR *A, TMatrix_DCSR *LD)
{
	return ERROR_NO_ERROR;
}

