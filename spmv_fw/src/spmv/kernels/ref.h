#ifndef _REF_H_
#define _REF_H_

#include <stdio.h>
#include <stdlib.h>
#include "spmv.h"

typedef struct {
    int hash;
    TMatrix_CSR matr;
} ref_desc_t;

#endif
