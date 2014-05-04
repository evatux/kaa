#include <stdio.h>
#include <stdlib.h>
#include "core.h"

int main(int argc, char **argv)
{
    int rows = (argc>1) ? atoi(argv[1]) : 4;
    int cols = (argc>2) ? atoi(argv[2]) : 6;

    int perm_rows[rows];
    int perm_cols[cols];

    /*
    perm_rows[rows-1] = rows-1;
    for (int i = 0; i < rows/2; ++i) {
        perm_rows[2*i] = 2*i + 1;
        perm_rows[2*i + 1] = 2*i;
    }

    perm_cols[cols-1] = cols-1;
    for (int i = 0; i < cols/2; ++i) {
        perm_cols[2*i] = 2*i + 1;
        perm_cols[2*i + 1] = 2*i;
    }
    */

    for (int i = 0; i < rows; ++i) perm_rows[i] = rows-1-i;
    for (int i = 0; i < cols; ++i) perm_cols[i] = cols-1-i;

    TMatrix_SMP matr_smp;
    matrix_smp_create(&matr_smp, rows, cols);
    for (int i = 0; i < rows*cols; ++i)
        matr_smp.val[i] = i;

    TMatrix_CSR matr1, matr2;
    matrix_smp2csr(&matr_smp, &matr1);
    matrix_perm_copy(&matr1, &matr2, perm_rows, perm_cols);

    matrix_show(&matr1, 1);
    matrix_show(&matr2, 1);

    return 0;
}
