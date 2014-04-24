#include "runit.h"

#include <stdlib.h>
#include <stdio.h>
#include "common.h"
#include "core.h"

static int is_pre_neps(TWGraph *gr, int v, int *dst, real threshold)
{
    *dst = -1;

    if (FABS(gr->wvert[v]) >= threshold)
        return 0;

    int flag, ci;
    int ng = gr->size;
    flag = 2;
    for (ci = gr->xadj[v]; ci < gr->xadj[v+1]; ci++)
    {
        if (gr->wedge[ci] != 0.) {
            if (flag==2) flag = 1;
            if (gr->adjncy[ci] < v) {
                flag = 0;
                break;
            } else {
                if (gr->adjncy[ci] < ng) ng = gr->adjncy[ci];
            }
        }
    }

    if (flag == 1) {
        *dst = ng;
        return 1;
    }

    return 0;
}

int graph_last_stage_reorder(TWGraph *gr, real threshold)
{
    if (threshold <= 0.) return ERROR_NO_ERROR;

    int *perm, *invp, *_invp, *neps_src, *neps_dst, *tt;
    int i, k, dst;
    int pre_neps_count = 0;
    int alpha, beta;
    int size = gr->size;
    int err = ERROR_NO_ERROR;

    perm = (int*)malloc(sizeof(int)*size);
    invp = (int*)malloc(sizeof(int)*size);
    _invp= (int*)malloc(sizeof(int)*size);
    neps_src = (int*)malloc(sizeof(int)*size);
    neps_dst = (int*)malloc(sizeof(int)*size);

    if ( perm == NULL || invp == NULL || _invp == NULL || neps_src == NULL || neps_dst == NULL ) {
        err = ERROR_MEMORY_ALLOCATION;
        goto finish;
    }

    for (i = 0; i < size; ++i)
    {
        if (is_pre_neps(gr, i, &dst, threshold))
        {
            neps_src[pre_neps_count] = i;
            neps_dst[pre_neps_count] = dst;
            pre_neps_count++;
        }
    }

    printf("pre_neps_count = %d\n", pre_neps_count);

    for (i = 0; i < size; ++i) invp[i] = i;

    for (i = 0; i < pre_neps_count; ++i)
    {
        alpha = invp[neps_src[i]];
        beta  = invp[neps_dst[i]];

        for (k = 0; k < size; ++k)
            _invp[k] = (alpha < invp[k] && invp[k] <= beta)?invp[k]-1:invp[k];
        _invp[neps_src[i]] = invp[neps_dst[i]];

        tt = invp; invp = _invp; _invp = tt;
   }

    for (i = 0; i < size; ++i) perm[invp[i]] = i;

    err = graph_reorder(gr, perm, invp);

finish:
    if (perm) free(perm);
    if (invp) free(invp);
    if (_invp) free(_invp);
    if (neps_src) free(neps_src);
    if (neps_dst) free(neps_dst);

    return err;
}
