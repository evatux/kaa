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
#define TOLERANCE    1e-5

int test_converter() {
    TMatrix_Simple m0, m1;
    TMatrix_DCSR   d1;
    int i, pass = 1;
    int size = MIN_TEST_DIM + rand() % (MAX_TEST_DIM - MIN_TEST_DIM);
    
    printf("TEST :: %10s :: size = %d\n", "converter", size);

    if ( (m0.val = (real*)malloc(sizeof(real)*size*size)) == NULL ) return ERROR_MEMORY_ALLOCATION;

    m0.size = size;
    for (i = 0; i < size*size; ++i) m0.val[i] = (1.0 / MAX_TEST_DIM)*(rand() % MAX_TEST_DIM);

    matrix_convert_simp2dcsr(&m0, &d1);
    matrix_convert_dcsr2simp(&d1, &m1);

    for (i = 0; i < size*size; ++i) if ( m0.val[i] != m1.val[i] ) { pass = 0; break; }

    matrix_simp_destroy(&m0);
    matrix_simp_destroy(&m1);
    matrix_destroy(&d1);

    if (!pass) printf("TEST :: %10s :: test failed (i = %d)", "converter", i);
    return (pass==1)?ERROR_NO_ERROR:19;
}

int test_solver() {
#define CLEANUP() \
    {                                   \
        free(X);                        \
        free(Y);                        \
        free(Z);                        \
        matrix_destroy(&LD);            \
        matrix_destroy(&A);             \
        matrix_simp_destroy(&simpM);    \
    }

    int i, j, p, non_zero_in_row;
    int neps;
    int err;
    int size = MIN_TEST_DIM + rand() % (MAX_TEST_DIM - MIN_TEST_DIM);
//  int size = 3;
    TMatrix_Simple simpM;
    TMatrix_DCSR   A, LD;
    real *X, *Y, *Z;
    real max_err, expected;

    printf("TEST :: %10s :: size = %d\n", "solver", size);

    simpM.size = size;
    if ( NULL == (simpM.val = (real*)malloc(sizeof(real)*size*size))) return ERROR_MEMORY_ALLOCATION;
    if ( NULL == (X         = (real*)malloc(sizeof(real)*size))     ) return ERROR_MEMORY_ALLOCATION;
    if ( NULL == (Y         = (real*)malloc(sizeof(real)*size))     ) return ERROR_MEMORY_ALLOCATION;
    if ( NULL == (Z         = (real*)malloc(sizeof(real)*size))     ) return ERROR_MEMORY_ALLOCATION;

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
        printf("TEST :: %10s :: neps = %d\n", "solver", neps);
        return 19;
    }

    for (i = 0; i < size; ++i) Z[i] = X[i] = SGN(rand() % 10 - 5) * rand() % ((int)(MAX_TEST_DIM * TEST_SPARSITY));

    err = solver(&LD, Z);
    if (err != ERROR_NO_ERROR) {
        printf("TEST :: %10s :: solver failed (err = %d)", "solver", err);
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
        printf("TEST :: %10s :: max_err = %.2e (expected: %.2e)\n", "solver",  max_err, expected);
        return 19;
    }
    return ERROR_NO_ERROR;
#undef CLEANUP
}

int test_lambda() {
    TMatrix_Simple SM;
    TMatrix_DCSR   A, LD;
    int i, neps, err;
    int size = 5; // MIN_TEST_DIM + rand() % (MAX_TEST_DIM - MIN_TEST_DIM);
    real M = -1, m = -1;
    
    printf("TEST :: %10s :: size = %d\n", "lambda", size);

    if ( (SM.val = (real*)malloc(sizeof(real)*size*size)) == NULL ) return ERROR_MEMORY_ALLOCATION;

    SM.size = size;
    for (i = 0; i < size*size; ++i) SM.val[i] = (i%size==i/size)?(i/size+1):0;
//  matrix_simp_show(&SM);

    err  = matrix_convert_simp2dcsr(&SM, &A);
    err |= cholesky_decomposition(&A, &LD, CHEPS_THRESHOLD, &neps);

    if (err != ERROR_NO_ERROR) {
        printf("TEST :: %10s :: preparation failed (err = %d)\n", "lambda", err);
        return err;
    }

//  matrix_show(&A, 1);
//  matrix_show(&LD, 1);

    err = lambda(&A, &LD, &m, &M);
    if (err != ERROR_NO_ERROR) {
        printf("TEST :: %10s :: test failed (err = %d)\n", "lambda", err);
        return err;
    }

    matrix_simp_destroy(&SM);
    matrix_destroy(&A);
    matrix_destroy(&LD);

    printf("TEST :: %10s :: l_min = %.2e, l_max = %.2e\n", "lambda", m, M);
    return ERROR_NO_ERROR;
}

int main(int argc, char** argv) {
    int err, i;
    srand( time(NULL) );

    if (argc >= 2) sscanf(argv[1], "%d", &i);
    switch(i) {
        case 1: goto l1; break;
        case 2: goto l2; break;
        case 3: goto l3; break;
    }

l1:
    err = test_converter();
    if (err != ERROR_NO_ERROR) printf("TEST :: %10s :: FAILED\n", "converter"), exit(1);

l2:
    err = test_solver();
    if (err != ERROR_NO_ERROR) printf("TEST :: %10s :: FAILED\n", "solver"), exit(1);

l3:
    err = test_lambda();
    if (err != ERROR_NO_ERROR) printf("TEST :: %10s :: FAILED\n", "lambda"), exit(1);

    printf("TEST PASSED\n");

    return 0;
}

