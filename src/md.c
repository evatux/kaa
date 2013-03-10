#include "md.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include "unit.h"

static void init_md(TMatrix_Simple* sd, int* degree, int *impacted)
{
    int size = sd->size;
    int i, j;
    for (i = 0; i < size; ++i)
    {
        degree[i]  = 0;
        impacted[i] = 0;
        for (j = 0; j < size; ++j)
        {      
            if (i != j && sd->val[i*size+j] != 0.) 
            {
                degree[i]++;
                sd->val[i*size+j] = 1.;
            }
        }
    }
}

// cp - current position
static int find_md(TMatrix_Simple* sd, int *degree, int *impacted, int at_least, int cp, real threshold)
{
    int size = sd->size;
    int i, j;
    int i0;             // answer

    int j0 = -1;
    real j0_val;        // --> min
    int j0_eps = -1;
    int j0_eps_imp;     // --> max

    int min_deg = -1;

    for (i = cp; i < size; ++i)
    {
        if ( ((min_deg >= degree[i]) || (min_deg == -1)) && (degree[i] > at_least) )
        {
            if (min_deg > degree[i] || min_deg == -1) {
                min_deg = degree[i];
                j0 = j0_eps = -1;
            }

            if ( FABS(sd->val[i*size+i]) > threshold )      // valueable vertex
            {
                if ( (j0 == -1) || (FABS(sd->val[i*size+i]) < j0_val) ) {
                    j0 = i;
                    j0_val = FABS(sd->val[i*size+i]);
                }
            }
            else                                            // epsilon vertex
            {
                if ( (j0_eps == -1) || ( impacted[i] > j0_eps_imp ) ) {
                    j0_eps = i;
                    j0_eps_imp = impacted[i];
                }
            }
        }
    }

    if ( min_deg == -1 ) return -2;                                 //  search is finished: at_least is too big

    if ( j0 != -1 ) return j0;                                      //  the prittiest answer

    if ( (j0_eps != -1) && (j0_eps_imp > 0) ) return j0_eps;        // not so good answer, but still oK

    if ( j0_eps != -1 ) {
        i0 = find_md(sd, degree, impacted, min_deg, cp, threshold); // search in level-up
        if ( i0 != -2 ) return i0;                                  // something found

        if ( (i0 == -2) && (at_least == -1) ) return j0_eps;        // if everything is bad
        return -2;                                                  // return low-level answer
    }
}

static void elements_switch(TMatrix_Simple *sd, int *degree, int *impacted, int e1, int e2)
{
    int size = sd->size;
    int i;

    if (e1 == e2) {
#ifdef _DEBUG_LEVEL_MD
        fprintf(stderr, "[md: elements_switch] e1 == e2 (%d)\n", e1);
#endif
        return;
    }

    if (degree)   SWAP(int, degree[e1],   degree[e2]);
    if (impacted) SWAP(int, impacted[e1], impacted[e2]);

    for (i = 0; i < size; ++i)
        SWAP(real, sd->val[e1*size+i], sd->val[e2*size+i]);

    for (i = 0; i < size; ++i)
        SWAP(real, sd->val[i*size+e1], sd->val[i*size+e2]);
}

static void elements_update(TMatrix_Simple *sd, int *degree, int *impacted, int el, real threshold)
{
    int size = sd->size;
    int i, j, jj;

    for (j = 0; j < size; ++j)
    {
        if (sd->val[el*size + j] != 0.)
        {
            if (FABS(sd->val[j*size + j]) < threshold) impacted[j]++;
            for (jj = j+1; jj < size; ++jj)
            {
                if (sd->val[el*size + jj] != 0.)
                {
                    if (sd->val[j*size + jj] != 0.) {
                        sd->val[j*size + jj] = sd->val[jj*size + j] = 1.;
                        degree[j]++;
                        degree[jj]++;
                    }
                }
            }
        }
    }
}

int md(TMatrix_DCSR *matr, real threshold)
{
    int err = ERROR_NO_ERROR;;
    TMatrix_Simple SD;
    int size;
    int i0;
    int cp = 0;
    int *degree   = NULL;
    int *impacted = NULL;
    int *perm = NULL, *invp = NULL;

#ifdef _DEBUG_LEVEL_MD
        fprintf(stderr, "[md: md] md started");
#endif

    size = matr->size;
    err = matrix_convert_dcsr2simp(matr, &SD);
    if (err != ERROR_NO_ERROR) goto finish;

    degree   = (int*)malloc( sizeof(int)*size );
    impacted = (int*)malloc( sizeof(int)*size );
    perm     = (int*)malloc( sizeof(int)*size );
    invp     = (int*)malloc( sizeof(int)*size );
    if (degree == NULL || impacted == NULL || perm == NULL || invp == NULL) {
        err = ERROR_MEMORY_ALLOCATION;
        goto finish;
    }

    init_md(&SD, degree, impacted);

#ifdef _DEBUG_LEVEL_MD
        fprintf(stderr, "[md: md] stage 1 started");
#endif
    while ( (i0 = find_md(&SD, degree, impacted, -1, cp, threshold) ) != -2 )
    {
        elements_switch(&SD, degree, impacted, i0, cp);
        elements_update(&SD, degree, impacted, cp, threshold);
        perm[cp] = i0;
        invp[i0] = cp;
        cp++;
    }

#ifdef _DEBUG_LEVEL_MD
    fprintf(stderr, "[md: md] stage 2 started\n");
#endif

    matrix_simp_destroy(&SD);
    err = matrix_convert_dcsr2simp(matr, &SD);
    if (err != ERROR_NO_ERROR) {
        fprintf(stderr, "[md: md] stage 2: matrix_convert_dcsr2simp failed\n");
        goto finish;
    }

    for (cp = 0; cp < size; cp++)
        elements_switch(&SD, NULL, NULL, cp, perm[cp]);

    matrix_destroy(matr);
    err = matrix_convert_simp2dcsr(&SD, matr);
    if (err != ERROR_NO_ERROR)
        fprintf(stderr, "[md: md] stage 2: matrix_convert_simp2dcsr failed\n");

#ifdef _DEBUG_LEVEL_MD
        fprintf(stderr, "[md: md] md finished");
#endif
finish:
    matrix_simp_destroy(&SD);
    if (degree)   free(degree);
    if (impacted) free(impacted);
    if (perm)     free(perm);
    if (invp)     free(invp);

    return err;
}

int md_mod (TMatrix_DCSR *matr)
{
    return md(matr, MD_EPS_THRESHOLD);
}

int md_orig(TMatrix_DCSR *matr)
{
    return md(matr, 0);
}
