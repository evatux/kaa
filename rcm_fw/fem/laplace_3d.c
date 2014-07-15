#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "common.h"
#include "core.h"

#define T_NAT 1
#define T_ND  2
#define T_QR  3

#define DIM 3

#define UNSET_VAL   (-1)
#define IS_UNSET(x) (x == UNSET_VAL)

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
    if (fabs(x - 0.3) < 0.2 || fabs(y - 0.6) < 0.2 || fabs(z - 0.2) < 0.1) return 100.;

    return 1.;
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

#define SET_IF_UNSET(x,y) \
    do { if (IS_UNSET(x)) (x) = (y); } while(0)

#define DEC(x) (--(*x))

void create_index_nd(int *ind, int n, int *last, int b[DIM], int e[DIM])
{
    typedef int (*ind_t)[n][n];
    ind_t X = (ind_t)ind;

    int c0, c1, c2, cc;

    for (cc = 0; cc < DIM; ++cc)
        if (b[cc] > e[cc]) return;

    int m[DIM], i[DIM];
    for (int c = 0; c < DIM; ++c)
        m[c] = (b[c] + e[c])/2;

    // dim 0
    X[m[0]][m[1]][m[2]] = DEC(last);

    // dim 1
    for (c0 = 0; c0 < DIM; ++c0)
    {
        for (cc = 0; cc < DIM; ++cc)
            if (cc != c0) i[cc] = m[cc];

        for (i[c0] = b[c0]; i[c0] < e[c0]; ++i[c0])
            SET_IF_UNSET(X[i[0]][i[1]][i[2]], DEC(last));
    }

    // dim 2
    for (c0 = 0; c0 < DIM; ++c0)
        for (c1 = c0+1; c1 < DIM; ++c1)
        {
            for (cc = 0; cc < DIM; ++cc)
                if (cc != c0 && cc != c1) i[cc] = m[cc];

            for (i[c0] = b[c0]; i[c0] < e[c0]; ++i[c0])
                for (i[c1] = b[c1]; i[c1] < e[c1]; ++i[c1])
                    SET_IF_UNSET(X[i[0]][i[1]][i[2]], DEC(last));
        }

    // the other :)
    int *bm[2] = { b, m };
    int *me[2] = { m, e };
    for (c0 = 0; c0 < 2; ++c0)
        for (c1 = 0; c1 < 2; ++c1)
            for (c2 = 0; c2 < 2; ++c2)
            {
                int bn[DIM] = {
                    bm[c0][0] + c0,
                    bm[c1][1] + c1,
                    bm[c2][2] + c2 };
                int en[DIM] = {
                    me[c0][0] - (1-c0),
                    me[c1][1] - (1-c1),
                    me[c2][2] - (1-c2) };
                create_index_nd(ind, n, last, bn, en);
            }
}

void create_index(int type, int *ind, int n, int last, int b[DIM], int e[DIM])
{
    switch(type) {
        case T_NAT: break;
        case T_ND: create_index_nd(ind, n, &last, b, e); break;
//        case T_QR: create_index_qr(ind, n, &last, b, e); break;
        default: fprintf(stderr, "wrong type\n"); exit(-2);
    }
}

int kij2e_nat(int n, int k, int i, int j)
{
    if (i <= 0 || i >= n) return -1;
    if (j <= 0 || j >= n) return -1;
    if (k <= 0 || k >= n) return -1;
    return ( (k-1) * (n-1) + (i-1) ) * (n-1) + (j-1);
}

int kij2e_ind(int *ind, int n, int k, int i, int j)
{
    typedef int (*ind_t)[n][n];
    ind_t X = (ind_t)ind;
    return X[k][i][j];
}

int kij2e(int type, int *ind, int n, int k, int i, int j)
{
    if (i <= 0 || i >= n) return -1;
    if (j <= 0 || j >= n) return -1;
    if (k <= 0 || k >= n) return -1;
    switch(type) {
        case T_NAT: return kij2e_nat(n, k, i, j);
        case T_ND:
        case T_QR:  return kij2e_ind(ind, n-1, k-1, i-1, j-1);
        default: fprintf(stderr, "wrong type\n"); exit(-2);
    }
}

static int get_n(const char *task, int *type)
{
    *type = T_NAT;
    if (task[0] >= '0' && task[0] <= '9') return atoi(task);
    switch (task[0]) {
        case 'n': *type = T_ND; break;
        case 'q': *type = T_QR; break;
        default:
                  fprintf(stderr, "incorrect type\n");
                  fprintf(stderr, "should be either 'n', or 'q' or nothing");
                  exit(-2);
    }
    return atoi(task+1);
}

int main(int argc, char **argv)
{
    int type;
    int n  = (argc > 1) ? get_n(argv[1], &type) : 10;
    const char* matr_fname = (argc > 2) ? argv[2] : "l3_default.csr";
    const char* pict_fname = (argc > 3) ? argv[3] : NULL;
    int *ind;

    real h  = 1. / n;
    int size = (n - 1) * (n - 1) * (n - 1);
    int nonz = 6 * size;

    int err;

    if (type != T_NAT) {
        ind = malloc(sizeof(int)*size);
        for (int i = 0; i < size; ++i) ind[i] = UNSET_VAL;

        int b[3] = { 0,   0,   0   };
        int e[3] = { n-1, n-1, n-1 };
        create_index(type, ind, n-1, size, b, e);
#if 0
        for(int k = 0; k < n-1; ++k, printf("\n"))
            for(int i = 0; i < n-1; ++i, printf("\n"))
                for(int j = 0; j < n-1; ++j)
                    printf("%-5d", ind[(k*(n-1)+i)*(n-1)+j]);
        exit(2);
#endif
    }

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
                    kij2e(type, ind, n, k - 1, i - 1, j - 1), // brl
                    kij2e(type, ind, n, k - 1, i - 1, j - 0), // brr
                    kij2e(type, ind, n, k - 1, i - 0, j - 1), // bfl
                    kij2e(type, ind, n, k - 1, i - 0, j - 0), // bfr
                    kij2e(type, ind, n, k - 0, i - 1, j - 1), // trl
                    kij2e(type, ind, n, k - 0, i - 1, j - 0), // trr
                    kij2e(type, ind, n, k - 0, i - 0, j - 1), // tfl
                    kij2e(type, ind, n, k - 0, i - 0, j - 0), // tfr
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
