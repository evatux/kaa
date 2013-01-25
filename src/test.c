#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "core.h"
#include "rcm.h"
#include "cholesky.h"
#include "solver.h"

#define MIN_TEST_DIM 50
#define MAX_TEST_DIM 1000
#define TEST_SPARSITY 0.1
#define TEST_CHEPS   1e-3
#define TOLERANCE	 1e-5

/*
int test_convert() {
	TMatrix_Simple simpM;
	TMatrix_DCSR   d1, d2;
	int i;

	matrix_load(&d1, config->matr_in_file);
	matrix_convert_dcsr2simp(&d1, &simpM);
	matrix_convert_simp2dcsr(&simpM, &d2);

#define CMI(X,I) do { if (d1.X != d2.X) { printf("test failed: (%d) %d != %d\n", I, d1.X, d2.X); exit(2); } } while(0)
#define CMF(X,I) do { if (d1.X != d2.X) { printf("test failed: (%d) %f != %f\n", I, d1.X, d2.X); exit(2); } } while(0)

	printf("\nsize: "); CMI(size,0);
	printf("\nnonz: "); CMI(nonz,0);
	printf("\ndiag: "); for (i = 0; i < d1.size; ++i) CMF(diag[i],i);
	printf("\nval:  "); for (i = 0; i < d1.nonz; ++i) CMF( val[i],i);
	printf("\ncol:  "); for (i = 0; i < d1.nonz; ++i) CMI(col_ind[i],i);
	printf("\nrow:  "); for (i = 0; i < d1.size; ++i) CMI(row_ptr[i],i);

	printf("test passed\n");
	exit(0);
}
*/
int test_solver() {
#define CLEANUP() \
	{									\
		free(X);						\
		free(Y);						\
		free(Z);						\
		matrix_destroy(&LD);			\
		matrix_destroy(&A);				\
		matrix_simp_destroy(&simpM);	\
	}

	int i, j, p, non_zero_in_row;
	int neps;
	int err;
	int size = MIN_TEST_DIM + rand() % (MAX_TEST_DIM - MIN_TEST_DIM);
//	int size = 3;
	TMatrix_Simple simpM;
	TMatrix_DCSR   A, LD;
	real *X, *Y, *Z;
	real max_err, expected;

	printf("TEST :: solver :: size = %d\n", size);

	simpM.size = size;
	if ( NULL == (simpM.val = (real*)malloc(sizeof(real)*size*size))) return ERROR_MEMORY_ALLOCATION;
	if ( NULL == (X 		= (real*)malloc(sizeof(real)*size)) 	) return ERROR_MEMORY_ALLOCATION;
	if ( NULL == (Y			= (real*)malloc(sizeof(real)*size)) 	) return ERROR_MEMORY_ALLOCATION;
	if ( NULL == (Z			= (real*)malloc(sizeof(real)*size)) 	) return ERROR_MEMORY_ALLOCATION;

/*
	simpM.val[0] = 1.; simpM.val[1] = 1.;simpM.val[2] = 0.;
	simpM.val[3] = 1.; simpM.val[4] = 2.;simpM.val[5] = 0.;
	simpM.val[6] = 0.; simpM.val[7] = 0.;simpM.val[8] = 1.;
*/

	// make SPD matrix
	for (i = 0; i < size*size; ++i) simpM.val[i] = 0;
	for (i = 0; i < size; ++i) simpM.val[i*size+i] = (size + 1) + rand() % (MAX_TEST_DIM - MIN_TEST_DIM);
	for (i = 1; i < size; ++i) {
		if (!(rand()%10)) continue; // lots of zero rows
		non_zero_in_row = 1 + rand() % ((int)(size * TEST_SPARSITY));
		for (p = 0; p < non_zero_in_row; ++p) {
			j = rand() % (i);
			simpM.val[j*size+i] = simpM.val[i*size+j] = SGN(rand() % 10 - 5) * (rand() % ((int)(size * TEST_SPARSITY)));
		}
	}

	matrix_convert_simp2dcsr(&simpM, &A);

	cholesky_decomposition(&A, &LD, TEST_CHEPS, &neps);
	if ( neps != 0 ) {
		CLEANUP();
		printf("TEST :: solver :: neps = %d\n", neps);
		return 19;
	}

	for (i = 0; i < size; ++i) Z[i] = X[i] = SGN(rand() % 10 - 5) * rand() % ((int)(MAX_TEST_DIM * TEST_SPARSITY));

	err = solver(&LD, Z);
	if (err != ERROR_NO_ERROR) {
		printf("TEST :: solver :: solver failed");
		CLEANUP();
		return 19;
	}

	matrix_vector_mult(&A, Z, Y);
/*
	printf("X\tY\tZ\t\tX ?= Y=AZ \n");
	for (i = 0; i < size; ++i) printf("%-8.2f%-8.2f%-8.2f\n", X[i], Y[i], Z[i]);
*/
	max_err = 0.; expected = X[0];
	for (i = 0; i < size; ++i) if (FABS(X[i] - Y[i]) > max_err) max_err = FABS(X[i] - Y[i]), expected = X[i];

	CLEANUP();

	if (max_err/(FABS(expected)+TOLERANCE) > TOLERANCE) {
		printf("TEST :: solver :: max_err = %.2e (expected: %.2e)\n", max_err, expected);
		return 19;
	}
	return ERROR_NO_ERROR;
#undef CLEANUP
}

int main(int argc, char** argv) {
	int err;
	srand( time(NULL) );

	err = test_solver();
	if (err != ERROR_NO_ERROR) printf("TEST :: solver :: FAILED\n"), exit(1);

	printf("TEST PASSED\n");

	return 0;
}
