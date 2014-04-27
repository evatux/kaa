#include <stdio.h>
#include <stdlib.h>
#include "core.h"
#include "spmv.h"

const int kernel_num = sizeof(kernel_list)/sizeof(kernel_list[0]);

double gflop(TMatrix_CSR *matr)
{
    return 1e-9*(2*matr->nonz - matr->rows);
}
