#ifndef _GBLOCK_OMP_H_
#define _GBLOCK_OMP_H_

#include <stdio.h>
#include <stdlib.h>
#include "spmv.h"

typedef struct {
    int hash;
    TMatrix_CSR matr;
} gblock_desc_t;

#endif
