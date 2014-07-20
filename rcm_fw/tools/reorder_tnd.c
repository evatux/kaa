#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "common.h"
#include "core.h"
#include "tnd.h"

#define ALG_ORIG 0
#define ALG_MULT 1

#define SAFE(f) \
    do { int err = f; \
        if (err != ERROR_NO_ERROR) \
        PRINT_ERROR_MESSAGE_AND_EXIT(err); \
    } while(0)

int in_row(TMatrix_DCSR *A, int row, int col)
{
    if (row == col && A->diag[row] != 0.) return 1;
    for (int i = A->row_ptr[row]; i < A->row_ptr[row+1]; ++i)
    {
        int cur_col = A->col_ind[i];
        if (cur_col <  col) continue;
        if (cur_col >  col) break;
        if (cur_col == col) return 1;
    }
    return 0;
}

int compare(const void *a, const void *b) { return (*(int*)a - *(int*)b); }

void do_orig(TMatrix_DCSR *A, int **_xadj, int **adjncy)
{
    int size = A->size;
    int total = 0;
    int *a;
    int *xadj = (int*)malloc(sizeof(int)*(size+1));
    if (!xadj)
        PRINT_ERROR_MESSAGE_AND_EXIT(ERROR_MEMORY_ALLOCATION);
    *_xadj = xadj;

    for (int step = 0; step <= 1; ++step)
    {
        if (step == 0) xadj[0] = total;
        if (step == 1)
        {
            a = (int*)malloc(sizeof(int)*total);
            if (!a) PRINT_ERROR_MESSAGE_AND_EXIT(ERROR_MEMORY_ALLOCATION);
            *adjncy = a;
        }

        for (int i = 0; i < size; ++i)
        {
            int cur = A->row_ptr[i+1] - A->row_ptr[i];

            // copy cur elements
            if (step == 1)
                for (int j = A->row_ptr[i]; j < A->row_ptr[i+1]; ++j)
                    a[xadj[i]+(j-A->row_ptr[i])] = A->col_ind[j];

            for (int j = 0; j < size; ++j)
            {
                if (!in_row(A,i,j) && in_row(A,j,i))
                {
                    if (step == 1) a[xadj[i] + cur] = j;
                    cur++;
                }
            }

            if (step == 0)
            {
                total += cur;
                xadj[i+1] = total;
            }
        }
    }

    for (int i = 0; i < size; ++i)
        qsort(&a[xadj[i]], xadj[i+1] - xadj[i], sizeof(int), compare);
}

void do_mult(TMatrix_DCSR *A, int **_xadj, int **adjncy)
{
    int size = A->size;
    int total = 0;
    int *a;
    int *xadj = (int*)malloc(sizeof(int)*(size+1));
    if (!xadj)
        PRINT_ERROR_MESSAGE_AND_EXIT(ERROR_MEMORY_ALLOCATION);
    *_xadj = xadj;
    xadj[0] = 0;

    for (int step = 0; step <= 1; ++step)
    {
        if (step == 1)
        {
            a = (int*)malloc(sizeof(int)*total);
            if (!a) PRINT_ERROR_MESSAGE_AND_EXIT(ERROR_MEMORY_ALLOCATION);
            *adjncy = a;
        }

        for (int i = 0; i < size; ++i)
        {
            int cur = 0;
            for (int j = 0; j < size; ++j)
            {
                if (i == j) continue;

                int is_there = 0;
                for (int ci = A->row_ptr[i]; ci < A->row_ptr[i+1]; ++ci)
                {
                    int cur_col = A->col_ind[ci];
                    if (in_row(A,j,cur_col)) {
                        is_there = 1;
                        break;
                    }
                }

                if (!is_there && in_row(A,j,i)) is_there = 1;

                if (is_there) {
                    if (step == 1) a[xadj[i] + cur] = j;
                    cur++;
                }
            }

            if (step == 0)
            {
                total += cur;
                xadj[i+1] = total;
            }
        }
    }
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "usage: %s infile outfile {o|m} [res_pic.png] [intermid_pic.png\n");
        return 2;
    }

    int err;
    int is_in_fmc  = (strstr(argv[1], ".mtx") != NULL);
    int is_out_fmc = (strstr(argv[2], ".mtx") != NULL);
    int type = (argv[3][0] == 'o') ? ALG_ORIG : ALG_MULT;

    TMatrix_DCSR _A,  *A  = &_A;
    TMatrix_DCSR _B,  *B  = &_B;
    TWGraph      _gr, *gr = &_gr;
    if (is_in_fmc) err = matrix_load_fmc(A, argv[1]);
    else err = matrix_load(A, argv[1]);
    if (err != ERROR_NO_ERROR)
        PRINT_ERROR_MESSAGE_AND_EXIT(err);

    int size = A->size;
    int *xadj, *adjncy;
    int *perm, *invp;

    if (!perm || !invp)
        PRINT_ERROR_MESSAGE_AND_EXIT(ERROR_MEMORY_ALLOCATION);

    switch (type) {
        case ALG_ORIG: do_orig(A, &xadj, &adjncy); break;
        case ALG_MULT: do_mult(A, &xadj, &adjncy); break;
    }

    TWGraph igr = { NULL, adjncy, xadj, NULL, size, xadj[size] };

    tnd_perm(&igr, &perm, &invp);

    SAFE(build_graph(gr, A));
    SAFE(graph_reorder(gr, perm, invp));
    SAFE(build_matrix(gr, B, 1));

    if (argc > 4) matrix_portrait(B, argv[4], 0, 0, NULL);
    if (is_out_fmc) matrix_save_fmc(B, argv[2]);
    else matrix_save(B, argv[2]);

    free(xadj);
    free(adjncy);
    if (argc > 5) {
        switch (type) {
            case ALG_ORIG: do_orig(B, &xadj, &adjncy); break;
            case ALG_MULT: do_mult(B, &xadj, &adjncy); break;
        }
        igr = (TWGraph){ NULL, adjncy, xadj, NULL, size, 0 };
        graph_portrait(&igr, argv[5]);
    }

    return 0;
}
