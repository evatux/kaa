#include "rcm.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include "unit.h"

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

static int find_periphery(TWGraph *gr, int *invp)
{
    int err = ERROR_NO_ERROR;
    int start_vertex = perm_get_free_index(gr->size, invp);
    if (start_vertex == -1) return -1;
    err = find_periphery_in_subgraph(gr, &start_vertex, NULL);
    if (err != ERROR_NO_ERROR) {
        fprintf(stderr, "{find_periphery}: find_periphery_in_subgraph exit with error code: %d\n", err);
        return -1;
    }
    return start_vertex;
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

    //  -2   -- not yet seen
    //  -1   -- already in queue
    //  0..n -- already numbered 

    err = perm_init(gr->size, _perm, _invp);
    if ( err != ERROR_NO_ERROR ) return err;

    perm = *_perm;
    invp = *_invp;

    int v, g, i, gg, gi, s = 0;
    int include_flag;
    TQueue queue;
    queue_init(&queue);

    while ( (root = find_periphery(gr, invp)) != -1 ) {

#ifdef _DEBUG_LEVEL_1
        printf("[debug(1)]{find_permutation}: root = %d\n", root);
#endif

        // actually it is simple bs
        invp[root] = -1;
        queue_push(&queue, root);

        while (queue_pop(&queue, &v)) {

            if (invp[v] != -1) continue;

            for (i = gr->xadj[v]; i < gr->xadj[v+1]; i++) {
                g = gr->adjncy[i];
                if ( ( invp[g] < 0 ) && ( FABS(gr->wvert[g]) < MIN2(threshold, FABS(gr->wvert[v])) ) ) {
                    // g - candidate for including in front of queue
#ifdef _ZRCM_SMART
                    include_flag = 1;
                    for (gi = gr->xadj[g]; gi < gr->xadj[g+1]; gi++) {
                        gg = gr->adjncy[gi];
                        if ( (invp[gg] < 0) && (gg != v) && ( FABS(gr->wvert[gg]) > MAX2(threshold, FABS(gr->wvert[g])) ) )
                        {
                            include_flag = 0;
                            break;
                        }
                    }

                    if (include_flag)
                    {
                        stack_push(&queue, v);
                        v = g;
                    }
                    break;
#else
                    stack_push(&queue, v);
                    v = g;
                    break;
#endif
                }
            }

#ifdef DIRECT_CM
            perm[s] = v;
            invp[v] = s++;                  // put connected vertices closer
#else
            perm[gr->size - 1 - s] = v;
            invp[v] = gr->size - 1 - s++;   // put connected vertices closer
#endif

            for (i = gr->xadj[v]; i < gr->xadj[v+1]; i++) {
                g = gr->adjncy[i];
                if ( invp[g] == -2 ) {
                    invp[g] = -1;
                    queue_push(&queue, g);  // put unindexed vertices to queue
                }
            }
        }
#ifdef _DEBUG_LEVEL_1
        printf("[debug]{find_permutation}: current permutation: \n\t");
        for (i = 0; i < gr->size; ++i) printf("%d ", perm[i]);
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

    return ERROR_NO_ERROR;
}

//  REORDERING: RCM. Reverse Cuthill-McKey

int rcm(TMatrix_DCSR *matr, real threshold)
{
    TWGraph gr;
    int *perm, *invp;
    int err = 0;

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{rcm}: graph_builder\n");
#endif
    err = build_graph(&gr, matr);
    if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("rcm: graph_builder failed", err);

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{rcm}: find_permutation\n");
#endif
    err = find_permutation(&gr,&perm, &invp, threshold);    // !!!REORDERING!!!
    if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("rcm: find_permutation failed", err);

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{rcm}: graph_reoder\n");
#endif
    err = graph_reorder(&gr, perm, invp);
    if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("rcm: graph_reorder failed", err);

    if (perm) free(perm);
    if (invp) free(invp);

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{rcm}: matrix_builder\n");
#endif
    err = build_matrix(&gr, matr, 0);
    if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("rcm: matrix_builder failed", err);

    graph_destroy(&gr);

    return ERROR_NO_ERROR;
}

int rcm_mod(TMatrix_DCSR *matr)
{
    return rcm(matr, EPS_THRESHOLD);
}

int rcm_orig(TMatrix_DCSR *matr)
{
    return rcm(matr, 0);
}

