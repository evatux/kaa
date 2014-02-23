#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "common.h"
#include "core.h"

#define brl 0
#define brr 1
#define bfl 2
#define bfr 3
#define trl 4
#define trr 5
#define tfl 6
#define tfr 7

// is_left, is_rear, is_bottom
#define is_l(x) (((x) & 1) == 0)
#define is_r(x) (((x) & 2) == 0)
#define is_b(x) (((x) & 4) == 0)

real get_coeff(real x, real y, real z)
{
    if (0 && fabs(x - 0.5) < 0.2 && fabs(y - 0.5) < 0.2 && fabs(z - 0.5) < 0.2) return 100.;

    return 1.;
}

int kij2e(int n, int k, int i, int j)
{
    if (i <= 0 || i >= n) return -1;
    if (j <= 0 || j >= n) return -1;
    if (k <= 0 || k >= n) return -1;
    return ( (k-1) * (n-1) + (i-1) ) * (n-1) + (j-1);
}

int loc_off(int c0, int c1)
{
    int diff = c1 - c0;

    switch (diff) {
        case -4: if (is_b(c1)) return 0; else break;  // bottom
        case -2: if (is_r(c1)) return 1; else break;  // rear
        case -1: if (is_l(c1)) return 2; else break;  // left
        case  1: if (is_l(c0)) return 3; else break;  // righ
        case  2: if (is_r(c0)) return 4; else break;  // forward
        case  4: if (is_b(c0)) return 5; else break;  // top
    }

    return -1;
}

int main(int argc, char **argv)
{
    int  n  = (argc > 1) ? atoi(argv[1]) : 10;
    const char* matr_fname = (argc > 2) ? argv[2] : "l3_default.csr";
    const char* pict_fname = (argc > 3) ? argv[3] : NULL;

    real h  = 1. / n;
    int size = (n - 1) * (n - 1) * (n - 1);
    int nonz = 6 * size;

    int err;

    TMatrix_DCSR _m;
    TMatrix_DCSR *m = &_m;
    err = matrix_create(m, size, nonz, 1);
    if (err) PRINT_ERROR_MESSAGE_AND_EXIT(err);

    for (int i = 1; i <= size; ++i)
        m->row_ptr[i] = 6 * i;

    const real esm[8][8] = {
        {  4, -1, -2,  0, -1,  0,  0,  0},
        { -1,  4,  0, -1,  0, -2,  0,  0},
        { -2,  0,  6, -2,  0,  0, -2,  0},
        {  0, -1, -2,  4,  0,  0,  0, -1},
        { -1,  0,  0,  0,  4, -2, -1,  0},
        {  0, -2,  0,  0, -2,  6,  0, -2},
        {  0,  0, -2,  0, -1,  0,  4, -1},
        {  0,  0,  0, -1,  0, -2, -1,  4}
    };

    for (int k = 1; k <= n; ++k)
    {
        for (int i = 1; i <= n; ++i)
        {
            for (int j = 1; j <= n; ++j)
            {
                int todo[8] = { 
                    kij2e(n, k - 1, i - 1, j - 1), // brl
                    kij2e(n, k - 1, i - 1, j - 0), // brr
                    kij2e(n, k - 1, i - 0, j - 1), // bfl
                    kij2e(n, k - 1, i - 0, j - 0), // bfr
                    kij2e(n, k - 0, i - 1, j - 1), // trl
                    kij2e(n, k - 0, i - 1, j - 0), // trr
                    kij2e(n, k - 0, i - 0, j - 1), // tfl
                    kij2e(n, k - 0, i - 0, j - 0), // tfr
                };

                for (int c0 = 0; c0 < 8; ++c0) {
                    int e = todo[c0];
                    if (e < 0) continue;
                    real coeff = h / 6. * get_coeff(j*h - h/2, i*h - h/2, k*h - h/2);
                    m->diag[e] += esm[c0][c0] * coeff;

                    for (int c1 = 0; c1 < 8; ++c1) {
                        int off = loc_off(c0, c1);
                        if (off < 0) continue;
                        m->col_ind[m->row_ptr[e] + off] =  todo[c1];
                        m->val    [m->row_ptr[e] + off] += esm[c0][c1] * coeff;
                    }
                }
            }
        }
    }

    TMatrix_DCSR _l3;
    TMatrix_DCSR *l3 = &_l3;
    matrix_copy_fix(m, l3);
    matrix_destroy(m);

    if (pict_fname) matrix_portrait(l3, pict_fname, 20., 0, NULL);
    matrix_save(l3, matr_fname);

    return 0;
}
