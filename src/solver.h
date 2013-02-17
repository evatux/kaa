#ifndef _SOLVER_H_
#define _SOLVER_H_

/***************************************************
* File: solver.h
* Set of matrix operation primitives:
*	o) w <- Bv, where B = A^-1, A - LDL^T
*	o) lambda {min, max} A
*	o) cond(A)
***************************************************/

#include "common.h"
#include "core.h"

#ifdef _TEST_DOUBLE
#define INV_EPS				1e-10
#define LAMBDA_MIN_DELTA	1e-4
#define LAMBDA_MAX_DELTA	1e-4
#define MAX_ITER			1000
#else
#define INV_EPS				1e-5
#define LAMBDA_MIN_DELTA	1e-2
#define LAMBDA_MAX_DELTA	1e-2
#define MAX_ITER			1000
#endif

int matrix_vector_mult(TMatrix_DCSR* /*A*/, real* /*X*/, real* /*Y*/);
int solver(TMatrix_DCSR* /*LD*/, real* /*Y*/);
int lambda(TMatrix_DCSR* /*A*/, TMatrix_DCSR* /*LD*/, real* /*l_min*/, real* /*l_max*/);

int matrix_conditinaly(TMatrix_DCSR* /*A*/, TMatrix_DCSR* /*LD*/);

#endif

