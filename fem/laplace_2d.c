#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "common.h"
#include "core.h"

#define bl 0
#define br 1
#define tl 2
#define tr 3

#define is_b(x) (((x) & 2) == 0)
#define is_l(x) (((x) & 1) == 0)

real get_coeff(real x, real y)
{
    if (fabs(x - 0.5) < 0.2 && fabs(y - 0.5) < 0.2) return 100.;

    return 1.;
}

int ij2k(int n, int i, int j)
{
    if (i <= 0 || i >= n) return -1;
    if (j <= 0 || j >= n) return -1;
    return (i - 1) * (n - 1) + (j - 1);
}

int loc_off(int ci, int cj)
{
    if (ci + cj == 3) return -1;

    int diff = cj - ci;

    switch (diff) {
        case -2: return 0;  // b
        case -1: return 1;  // l
        case  1: return 2;  // r
        case  2: return 3;  // t
    }

    return -1;
}

int main(int argc, char **argv)
{
    int  n  = (argc > 1) ? atoi(argv[1]) : 100;
    const char* matr_fname = (argc > 2) ? argv[2] : "l2_default.csr";
    const char* pict_fname = (argc > 3) ? argv[3] : NULL;

    real h  = 1. / n;
    int size = (n - 1) * (n - 1);
    int nonz = 4 * size;

    int err;

    TMatrix_DCSR _m;
    TMatrix_DCSR *m = &_m;
    err = matrix_create(m, size, nonz, 1);
    if (err) PRINT_ERROR_MESSAGE_AND_EXIT(err);

    for (int i = 1; i <= size; ++i)
        m->row_ptr[i] = 4 * i;

    for (int i = 1; i <= n; ++i)
    {
        for (int j = 1; j <= n; ++j)
        {
            int todo[4] = { 
                            ij2k(n, i - 1, j - 1), // bl
                            ij2k(n, i - 1, j - 0), // br
                            ij2k(n, i - 0, j - 1), // tl
                            ij2k(n, i - 0, j - 0), // tr
                          };

            for (int ci = 0; ci < 4; ++ci) {
                int k = todo[ci];
                if (k < 0) continue;
                real coeff = get_coeff(j*h - h/2, i*h - h/2);
                m->diag[k] += 1 * coeff;

                for (int cj = 0; cj < 4; ++cj) {
                    int off = loc_off(ci, cj);
                    if (off < 0) continue;
                    m->col_ind[m->row_ptr[k] + off] =  todo[cj];
                    m->val    [m->row_ptr[k] + off] -= 0.5 * coeff;
                }
            }
        }
    }

    TMatrix_DCSR _l2;
    TMatrix_DCSR *l2 = &_l2;
    matrix_copy_fix(m, l2);
    matrix_destroy(m);

    if (pict_fname) matrix_portrait(l2, pict_fname, 5., 0, NULL);
    matrix_save(l2, matr_fname);

    return 0;
}
