#include "wnd.h"
#include "nd.h"

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

static int weight_map(TWGraph *gr, real **_map)
{
    int size = gr->size;
    *_map = (real*)malloc(sizeof(real)*size);
    if (NULL == *_map) return ERROR_MEMORY_ALLOCATION;
    real *map = *_map;

    int g, v;
    int level;
    int next_level_g;
    for (g = 0; g < size; ++g)
    {
        TQueue queue;
        queue_init(&queue);
        queue_push(&queue, g);
        level = 0;
        next_level_g = 0;

        while (queue_pop(&queue, &v)) {
            if (v == next_level_g) {
                ++level;
                next_level

        }
    }
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

static int nd_for_subgraph(TWGraph *gr, int start_vertex, int *perm, int *invp, int *N, int *ind_ptr, int *level, real threshold)
{
    int size = gr->size;
    int root = find_periphery(gr, start_vertex, invp);
    int n, g, l;
    int err;

    if (root == -1) return ERROR_SMTH_WRONG;

    graph_level(gr, root, invp, ind_ptr, level, &n);
    if (n <= 3) {   //  S <-- C
        int cl;
            
        //  first include small elements
        for (cl = 0; cl < n; cl++)
        {
            for (l = ind_ptr[cl]; l < ind_ptr[cl+1]; l++)
            {
                g = level[l];
                if (FABS(gr->wvert[g]) >= threshold) continue;
                (*N)-=1;
                invp[g] = *N;
                perm[*N] = g;
            }
        }
        //  then the others
        for (cl = 0; cl < n; cl++)
        {
            for (l = ind_ptr[cl]; l < ind_ptr[cl+1]; l++)
            {
                g = level[l];
                if (FABS(gr->wvert[g]) < threshold) continue;
                (*N)-=1;
                invp[g] = *N;
                perm[*N] = g;
            }
        }
    } else {          //  S <-- L[n/2+1]
        int j;
        int include_flag, include_flag_eps;
        int g_eps;
        int ii, jj;

#if (ND_SMART_TYPE == 1)
        {
            j = n/2;

            for (l = ind_ptr[j]; l < ind_ptr[j+1]; l++) {
                include_flag = 0;
                g = level[l];

                if (FABS(gr->wvert[g]) < threshold) {
                    include_flag = 1;
                } else {
                    for (ii = gr->xadj[g]; ii < gr->xadj[g+1]; ii++)
                    {
                        for (jj = ind_ptr[j+1]; jj < ind_ptr[j+2]; jj++)
                        {
                            if (gr->adjncy[ii] == level[jj]) {
                                include_flag = 1;
                                break;
                            }
                        }
                        if (include_flag) break;
                    }
                }

                if (include_flag) {
                    (*N)-=1;
                    invp[g] = *N;
                    perm[*N] = g;
                }
            }
        }
#elif (ND_SMART_TYPE == 2)        
        {
            j = n/2 - ((n&1)?0:1); 
            for (l = ind_ptr[j]; l < ind_ptr[j+1]; l++) {
                include_flag = 0;
                g = level[l];

                if (FABS(gr->wvert[g]) < threshold) {
                    include_flag = 1;
                } else {
                    for (ii = gr->xadj[g]; ii < gr->xadj[g+1]; ii++)
                    {
                        for (jj = ind_ptr[j+1]; jj < ind_ptr[j+2]; jj++)
                        {
                            if (gr->adjncy[ii] == level[jj]) {
                                include_flag = 1;
                                break;
                            }
                        }
                        if (include_flag) break;
                    }
                }

                if (include_flag) {
                    (*N)-=1;
                    invp[g] = *N;
                    perm[*N] = g;
                }
            }

            for (l = ind_ptr[j+1]; l < ind_ptr[j+2]; l++) {
                include_flag = 0;
                g = level[l];

                if (FABS(gr->wvert[g]) < threshold) {
                    include_flag = 1;
                }

                if (include_flag) {
                    (*N)-=1;
                    invp[g] = *N;
                    perm[*N] = g;
                }
            }
        }
#endif

        int g1 = level[0];
        int g2 = level[ind_ptr[n] - 1];

#ifdef _DEBUG_LEVEL_1
        int __i;
        printf("[nd_1]{find_permutation}: current permutation: \n\t");
        for (__i = 0; __i < gr->size; ++__i) printf("%d ", perm[__i]);
        printf("\n");
#endif

        err = nd_for_subgraph(gr, g1, perm, invp, N, ind_ptr, level, threshold);
        if (err != ERROR_NO_ERROR) ERROR_MESSAGE("wnd: nd_for_subgraph_1 failed", err);

#ifdef _DEBUG_LEVEL_1
        printf("[nd_2]{find_permutation}: current permutation: \n\t");
        for (__i = 0; __i < gr->size; ++__i) printf("%d ", perm[__i]);
        printf("\n");
#endif

        err = nd_for_subgraph(gr, g2, perm, invp, N, ind_ptr, level, threshold);
        if (err != ERROR_NO_ERROR) ERROR_MESSAGE("wnd: nd_for_subgraph_2 failed", err);
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

        err = nd_for_subgraph(gr, root, perm, invp, &n, ind_ptr, level, threshold);
        if (err != ERROR_NO_ERROR) ERROR_MESSAGE("wnd: nd_for_subgraph failed", err);

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

static int wnd_for_subgraph(TWGraph *gr, int start_vertex, int *perm, int *invp, int *N, int *ind_ptr, int *level, real threshold)
{
    int n    = gr->size; //*N;
    int g, l;
    int err = ERROR_NO_ERROR;

    real min = 1.e100;
    real max = 0.;

    static int *visited = NULL;
    if (!visited) visited = (int*)malloc(sizeof(int)*n);
    for (g = 0; g < n; ++g) visited[g] = 0;

    TQueue queue, queue_clique;
    queue_init(&queue);
    queue_push(&queue, start_vertex);

    while (queue_pop(&queue, &g))
    {
        if (visited[g]) continue;
        visited[g] = 1;

        real aw = FABS(gr->wvert[g]);
        if (aw < min) min = aw;
        if (aw > max) max = aw;

        for (l = gr->xadj[g]; l < gr->xadj[g+1]; ++l) {
            int h = gr->adjncy[l];  
            if (!visited[h] && invp[h] == -2) queue_push(&queue, h);
        }
    }

    if ( max < WND_AM_THRESHOLD * min )
    {
        //  the difference between min and max weights is quite small
        //  so we can switch to normal ND algorithm

        free(visited); visited = NULL;
        int err = ERROR_NO_ERROR;
        err = nd_for_subgraph(gr, start_vertex, perm, invp, N, ind_ptr, level, 0.);
        return err;
    }

    // let's divide by weight value
    // cut by max value
    real mean = (max + min) / 2.;

    for (g = 0; g < n; ++g) visited[g] = -2;
    queue_init(&queue);
    queue_init(&queue_clique);
    queue_push(&queue, start_vertex);

    int max_clique = 0;
    int max_root   = -1;
    while (queue_pop(&queue, &g))
    {
        if (visited[g] >= 0) continue;
        if (FABS(gr->wvert[g]) > mean)
        {
            int cur_vnum = 0;
            int v;
            queue_push(&queue_clique, g);
            while (queue_pop(&queue_clique, &v))
            {
                cur_vnum++;
                visited[v] = g;
                for (l = gr->xadj[v]; l < gr->xadj[v+1]; ++l)
                {
                    int h = gr->adjncy[l];
                    if (visited[h] == -2 && invp[h] == -2)
                    {
                        if (FABS(gr->wvert[h]) > mean) queue_push(&queue_clique, h);
                        else queue_push(&queue, h);
                    }
                }
            }

            if (cur_vnum > max_clique)
            {
                max_clique = cur_vnum;
                max_root   = g;
            }
        } else /* FABS <= mean */
        {
            if (visited[g] == -1) continue;
            visited[g] = -1;

            for (l = gr->xadj[g]; l < gr->xadj[g+1]; ++l)
            {
                int h = gr->adjncy[l];
                if (visited[h] == -2 && invp[h] == -2)
                    queue_push(&queue, h);
            }
        }
    }

    if (max_root == -1)
    {
        fprintf(stderr, "[debug]{wnd_for_subgraph}: max clique isn't found\n");
        exit(100);
    }

    queue_init(&queue);
    queue_init(&queue_clique);
    for (g = 0; g < n; ++g) visited[g] = 0;
    queue_push(&queue_clique, max_root);

    int first_g = -2;

    while(queue_pop(&queue_clique, &g))
    {
        if (visited[g]) continue;

//#define WND_AM_OUT_RING_SEPARATOR
#ifdef  WND_AM_OUT_RING_SEPARATOR
        if (FABS(gr->wvert[g]) < mean)
        {
            if (invp[g] == -2)
            {
                (*N)-=1;
                invp[g] = *N;
                perm[*N] = g;
            }
            continue;            
        }
        for (l = gr->xadj[g]; l < gr->xadj[g+1]; ++l)
        {
            int h = gr->adjncy[l];
            if (!visited[h])
                queue_push(&queue_clique, h);
        }
#else
        int include_flag = 0;
        for (l = gr->xadj[g]; l < gr->xadj[g+1]; ++l)
        {
            int h = gr->adjncy[l];
            if (visited[h] || invp[h] != -2) continue;

            if (FABS(gr->wvert[h]) > mean)
                queue_push(&queue_clique, h);
            else {
                queue_push(&queue, h);
                include_flag = 1;
            }
        }

        if (include_flag)
        {
            (*N)-=1;
            invp[g] = *N;
            perm[*N] = g;
        } else
        {
            if (first_g == -2) first_g = g;
        }

        visited[g] = 1;
#endif
    }

    free(visited); visited = NULL;

    if (first_g != -2) wnd_for_subgraph(gr, first_g, perm, invp, N, ind_ptr, level, threshold);
    while (queue_pop(&queue, &g))
    {
        if (invp[g] != -2) continue;
        wnd_for_subgraph(gr, g, perm, invp, N, ind_ptr, level, threshold);
    }

    return err;
}

// find permuation using arithmetical mean dividing by 2 (greater than am, and less than)
static int find_permutation_am(TWGraph *gr, int **_perm, int **_invp, real threshold)
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
        printf("[debug(1)]{find_permutation_am}: root = %d\n", root);
#endif

        err = wnd_for_subgraph(gr, root, perm, invp, &n, ind_ptr, level, threshold);
        if (err != ERROR_NO_ERROR) ERROR_MESSAGE("wnd: wnd_for_subgraph failed", err);

#ifdef _DEBUG_LEVEL_1
        int __i;
        printf("[debug]{find_permutation_am}: current permutation: \n\t");
        for (__i = 0; __i < gr->size; ++__i) printf("%d ", perm[__i]);
        printf("\n");
#endif
    }

#ifdef _DEBUG_LEVEL_1
    int __i;
    printf("[debug(1)]{find_permutation_am}:");
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

//  REORDERING: WND. Weight Nested Dissection Method

int wnd(TMatrix_DCSR *matr, real threshold)
{
    TWGraph gr;
    int *perm, *invp;
    int err = 0;

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{wnd}: graph_builder\n");
#endif
    err = build_graph(&gr, matr);
    if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("wnd: graph_builder failed", err);

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{wnd}: find_permutation\n");
#endif
    err = find_permutation_am(&gr,&perm, &invp, threshold);    // !!!REORDERING!!!
    if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("wnd: find_permutation failed", err);

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{wnd}: graph_reoder\n");
#endif
    err = graph_reorder(&gr, perm, invp);
    if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("wnd: graph_reorder failed", err);

    if (perm) free(perm);
    if (invp) free(invp);

    err = graph_last_stage_reorder(&gr, threshold);
    if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("wnd: graph_last_stage_reorder failed", err);

#ifdef _DEBUG_LEVEL_0
printf("[debug(0)]{wnd}: matrix_builder\n");
#endif
    err = build_matrix(&gr, matr, 0);
    if ( err != ERROR_NO_ERROR ) ERROR_MESSAGE("wnd: matrix_builder failed", err);

    graph_destroy(&gr);

    return ERROR_NO_ERROR;
}

int wnd_mod(TMatrix_DCSR *matr)
{
    return wnd(matr, WND_EPS_THRESHOLD);
}

int wnd_orig(TMatrix_DCSR *matr)
{
    return wnd(matr, 0);
}

