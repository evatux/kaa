#ifndef _CHOLESKY_H_
#define _CHOLESKY_H_

/***************************************************
* File: colesky.h
* Realization Cholesky decomposition
*
* NOTE: requires sizeof(real)*N**2 memory
***************************************************/

#include "common.h"
#include "core.h"

#define CHEPS_THRESHOLD		1e-3

//	Cholesky decomposition interface
int cholesky_decomposition(TMatrix_DCSR* /*A*/, TMatrix_DCSR* /*LD*/, const real /*cheps*/, int* /*neps*/);

#endif
