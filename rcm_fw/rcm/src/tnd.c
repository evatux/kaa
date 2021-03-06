#include "tnd.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include "unit.h"
#include "runit.h"

#define ERROR_MESSAGE(c,e) do { printf("%s (err_code: %d)\n", c, e); return e; } while(0)

static int perm_init(int size, int **_perm, int **_invp)
{
    int i;
    int *perm = (int*)malloc(sizeof(int)*size);
    int *invp = (int*)malloc(sizeof(int)*size);

    if ( NULL == perm || NULL == invp ) {
        if (perm) free(perm);
        if (invp) free(invp);
        *_perm = *_invp = NULL;
        return ERROR_MEMORY_ALLOCATION;
    }

    for (i = 0; i < size; ++i) perm[i] = invp[i] = -2;

    *_perm = perm;
    *_invp = invp;

    return ERROR_NO_ERROR;
}

static int perm_get_free_index(int size, int *invp)
{
    int i;

    for (i = 0; i < size; ++i) if ( invp[i] == -2 ) return i;
                          else if ( invp[i] < 0 ) printf("[debug]: {perm_get_free_index}: invp[%d] = %d\n", i, invp[i]), exit(2);
#ifdef _DEBUG_LEVEL_0
                          else if ( invp[i] < 0 ) printf("[debug]: {perm_get_free_index}: invp[%d] = %d\n", i, invp[i]), exit(3);
#endif

    return -1;
}

static int find_periphery(TWGraph *gr, int start_vertex, int* excepted)
{
    int err = ERROR_NO_ERROR;
    int vertex = start_vertex;
    err = find_periphery_in_subgraph(gr, &vertex, excepted);
    if (err != ERROR_NO_ERROR) {
        fprintf(stderr, "{find_periphery}: find_periphery_in_subgraph exit with error code: %d\n", err);
        return -1;
    }
    return vertex;
}

static int tnd_for_subgraph(TWGraph *gr, int start_vertex, int *perm, int *invp, int *N, int *ind_ptr, int *level, real threshold)
{
    int size = gr->size;
    int root = find_periphery(gr, start_vertex, invp);
    int n, g, l;
    int err;

    if (root == -1) return ERROR_SMTH_WRONG;

    graph_level(gr, root, invp, ind_ptr, level, &n);
    if (n <= 3) {   //  S <-- C
        int cl;

        for (cl = n/2+(n&1); cl < n; cl++)
        {
            for (l = ind_ptr[cl]; l < ind_ptr[cl+1]; l++)
            {
                g = level[l];
                (*N) -= 1;
                invp[g]  = *N;
                perm[*N] =  g;
            }
        }
        for (cl = n/2+(n&1)-1; cl >= 0; cl--)
        {
            for (l = ind_ptr[cl]; l < ind_ptr[cl+1]; l++)
            {
                g = level[l];
                (*N) -= 1;
                invp[g]  = *N;
                perm[*N] =  g;
            }
        }
    } else
    {               //  S <-- L[n/2+1]
        int jj, j, j0;
        j0 = n/2 + (n&1);

        {
            int map[3] = { 0, -1, 1 };

            for (jj = 0; jj < 3; ++jj)
            {
                j = j0 + map[jj];
                if (j >= n || j < 0) continue;

                for (l = ind_ptr[j]; l < ind_ptr[j+1]; l++) {
                    g = level[l];
                    (*N) -= 1;
                    invp[g]  = *N;
                    perm[*N] =  g;
                }
            }
        }

        int g1 = level[0];
        int g2 = level[ind_ptr[n] - 1];

#ifdef _DEBUG_LEVEL_1
        int __i;
        printf("[tnd_1]{find_permutation}: current permutation: \n\t");
        for (__i = 0; __i < gr->size; ++__i) printf("%d ", invp[__i]);
        printf("\n");
#endif

        if (invp[g1] == -2)
        {
            err = tnd_for_subgraph(gr, g1, perm, invp, N, ind_ptr, level, threshold);
            if (err != ERROR_NO_ERROR) ERROR_MESSAGE("nd: nd_for_subgraph_1 failed", err);
        }

#ifdef _DEBUG_LEVEL_1
        printf("[tnd_2]{find_permutation}: current permutation: \n\t");
        for (__i = 0; __i < gr->size; ++__i) printf("%d ", invp[__i]);
        printf("\n");
#endif

        if (invp[g1] == -2)
        {
            err = tnd_for_subgraph(gr, g2, perm, invp, N, ind_ptr, level, threshold);
            if (err != ERROR_NO_ERROR) ERROR_MESSAGE("tnd: tnd_for_subgraph_2 failed", err);
        }
    }

    return ERROR_NO_ERROR;
}

static int find_permutation(TWGraph *gr, int **_perm, int **_invp, real threshold)
{
// Making permututation vector
//
// v-th vertex comes perm[v]
// e.g. perm[0] <-- root
//      invp[root] <-- 0

    int err = ERROR_NO_ERROR;
    int root;
    int *perm, *invp;

    int *level, *ind_ptr;
    int n = gr->size;

    level   = (int*)malloc(sizeof(int)*(gr->size));
    ind_ptr = (int*)malloc(sizeof(int)*(gr->size+1));

    if (level == NULL || ind_ptr == NULL) {
        if (level) free(level);
        if (ind_ptr) free(ind_ptr);
        return ERROR_MEMORY_ALLOCATION;
    }

    //  -2   -- not yet seen
    //  -1   -- already in queue
    //  0..n -- already numbered

    err = perm_init(gr->size, _perm, _invp);
    if ( err != ERROR_NO_ERROR ) return err;

    perm = *_perm;
    invp = *_invp;

    while ( (root = perm_get_free_index(gr->size, invp)) != -1 ) {
#ifdef _DEBUG_LEVEL_1
        printf("[debug(1)]{find_permutation}: root = %d\n", root);
#endif

        err = tnd_for_subgraph(gr, root, perm, invp, &n, ind_ptr, level, threshold);
        if (err != ERROR_NO_ERROR) ERROR_MESSAGE("tnd: tnd_for_subgraph failed", err);

#ifdef _DEBUG_LEVEL_1
        int __i;
        printf("[debug]{find_permutation}: current permutation: \n\t");
        for (__i = 0; __i < gr->size; ++__i) printf("%d ", perm[__i]);
        printf("\n");
#endif
    }

#ifdef _DEBUG_LEVEL_1
    int __i;
    printf("[debug(1)]{find_permutation}:");
    printf("\n   i:");
    for (__i = 0; __i < gr->size; __i++) printf("%4d ", __i);
    printf("\nperm:");
    for (__i = 0; __i < gr->size; __i++) printf("%4d ", perm[__i]);
    printf("\ninvp:");
    for (__i = 0; __i < gr->size; __i++) printf("%4d ", invp[__i]);
    printf("\n");
#endif

    if (level) free(level);
    if (ind_ptr) free(ind_ptr);

    return ERROR_NO_ERROR;
}

//  REORDERING: TND. Triple Nested Dissection Method

int tnd(TMatrix_DCSR *matr, real threshold)
{
    TWGraph gr;
    int *perm, *invp;
    int err = 0;

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{tnd}: graph_builder\n");
#endif
    err = build_graph(&gr, matr);
    if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("tnd: graph_builder failed", err);

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{tnd}: find_permutation\n");
#endif
    err = find_permutation(&gr,&perm, &invp, threshold);    // !!!REORDERING!!!
    if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("tnd: find_permutation failed", err);

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{tnd}: graph_reoder\n");
#endif
    err = graph_reorder(&gr, perm, invp);
    if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("tnd: graph_reorder failed", err);

    if (perm) free(perm);
    if (invp) free(invp);

    err = graph_last_stage_reorder(&gr, threshold);
    if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("tnd: graph_last_stage_reorder failed", err);

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{tnd}: matrix_builder\n");
#endif
    err = build_matrix(&gr, matr, 0);
    if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("tnd: matrix_builder failed", err);

    graph_destroy(&gr);

    return ERROR_NO_ERROR;
}

int tnd_perm(TWGraph *gr, int **_perm, int **_invp)
{
    int err = 0;

    err = find_permutation(gr, _perm, _invp, 0);
    if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("tnd: find_permutation failed", err);

    return ERROR_NO_ERROR;
}

int tnd_mod(TMatrix_DCSR *matr)
{
    return tnd(matr, TND_EPS_THRESHOLD);
}

int tnd_orig(TMatrix_DCSR *matr)
{
    return tnd(matr, 0);
}

