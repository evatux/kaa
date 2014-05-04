#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "core.h"
#include "spmv.h"

const int kernel_num = sizeof(kernel_list)/sizeof(kernel_list[0]);

const char* ker_get_opt(info_t info, const char c)
{
    static char s[MAX_INFO_SIZE];
    strcpy(s, (const char*)info);
    char *pch = strtok(s, ":");
    while (pch != NULL) {
        if (pch[0] == c) {
            strcpy(s, pch+1);
            return s;
        }
        pch = strtok(NULL, ":");
    }
    return NULL;
}

double gflop(TMatrix_CSR *matr)
{
    return 1e-9*(2*matr->nonz - matr->rows);
}
