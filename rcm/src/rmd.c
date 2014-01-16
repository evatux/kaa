#include "rmd.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include "unit.h"

#define ERROR_MESSAGE(c,e) do { printf("%s (err_code: %d)\n", c, e); return e; } while(0)

static int compare(const void *a, const void *b)
{
    return (*(int*)a - *(int*)b);
}

static int find_permutation(TWGraph *gr, int **_perm, int **_invp, real threshold)
{
// Making permutation vector
//
// v-th vertex comes perm[v] 
// e.g. perm[0] <-- root
//      invp[root] <-- 0

    int err = ERROR_NO_ERROR;
    int *perm, *invp;
    int i, size = gr->size;
    int n = (int)threshold;  // FIXME: of course it is not threshold, but level of depth :)
    int v, w;

    perm = (int*)malloc(sizeof(int)*size);
    invp = (int*)malloc(sizeof(int)*size);
    if ( !perm || !invp )
    {
        if (perm) free(perm);
        if (invp) free(invp);
        return ERROR_MEMORY_ALLOCATION;
    }
    *_perm = perm;
    *_invp = invp;

    // vdegr[0] = { ... } - auxilary ((n-1)-th depth) "degree"
    // vdegr[1] = { ... } - final    (    n-th depth) "degree"
    int *array = (int*)malloc(sizeof(int)*2*size);
    if ( !array )
    {
        free(perm);
        free(invp);
        return ERROR_MEMORY_ALLOCATION;
    }

    int *vdegr[2] = { &array[0], &array[size] };

    for (v = 0; v < size; ++v)
    {
        vdegr[0][v] = gr->xadj[v+1] - gr->xadj[v];
        vdegr[1][v] = vdegr[0][v];
    }

    // level of depth
    for (i = 1; i < n; ++i)
    {
        // foreach (v) 
        //   foreach(w from Adj(v) || w = v)
        //     degree`[v] += degree[w];
        for (v = 0; v < size; ++v)
        {
            for (w = gr->xadj[v]; w < gr->xadj[v+1]; ++w)
                vdegr[1][v] += vdegr[0][gr->adjncy[w]];
        }

        // updating level
        for (v = 0; v < size; ++v)
            vdegr[0][v] = vdegr[1][v];
    }

    // permutation is simple sorted vdegr[1]-array 
    // from smaller to larger value of "degree"
    
    // make an array ( (degr[v], v) )
    for (v = 0; v < size; ++v)  // this cycle has data dependencies
    {
        array[2*v+0] = vdegr[1][v];
        array[2*v+1] = v;
    }

    qsort(array, size, 2*sizeof(int), compare); // hack: sort pairs (degr[v], v) by first argument

    for (v = 0; v < size; ++v)
    {
        w = array[2*v+1]; // w --> v
        perm[w] = v;
        invp[v] = w;
    }

    free(array);

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

    return ERROR_NO_ERROR;
}

//  REORDERING: RMD. Reverse Minimum Degree

int rmd(TMatrix_DCSR *matr, real threshold)
{
    TWGraph gr;
    int *perm, *invp;
    int err = ERROR_NO_ERROR;

    if (threshold < 1) return ERROR_NO_ERROR;

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{rmd}: graph_builder\n");
#endif
    err = build_graph(&gr, matr);
    if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("rmd: graph_builder failed", err);

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{rmd}: find_permutation\n");
#endif
    err = find_permutation(&gr,&perm, &invp, threshold);    // !!!REORDERING!!!
    if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("rmd: find_permutation failed", err);

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{rmd}: graph_reoder\n");
#endif
    err = graph_reorder(&gr, perm, invp);
    if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("rmd: graph_reorder failed", err);

    if (perm) free(perm);
    if (invp) free(invp);

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{rmd}: matrix_builder\n");
#endif
    err = build_matrix(&gr, matr, 0);
    if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("rmd: matrix_builder failed", err);

    graph_destroy(&gr);

    return ERROR_NO_ERROR;
}

int rmd_mod(TMatrix_DCSR *matr)
{
    return rmd(matr, RMD_EPS_THRESHOLD);
}

int rmd_orig(TMatrix_DCSR *matr)
{
    return rmd(matr, 0);
}

