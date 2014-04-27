#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include "common.h"
#include "core.h"
#include "spmv.h"
#include "bench.h"

double timer()
{
    struct timeval tm;
    gettimeofday( &tm, NULL );
    return (double)tm.tv_sec + (double)tm.tv_usec / 1e6;
}

int perf(spmv_kernel_t *ker, info_t info, perf_t *p)
{
    TMatrix_CSR *matr = p->matr;
    int nnz     = p->matr->nonz * p->rhs;
    int batch   = p->batch ? p->batch :
        ((nnz > PERF_BIG_NNZ) ? PERF_BIG_BATCH : PERF_SMALL_BATCH);
    int rounds  = p->rounds ? p->rounds :
        ((nnz > PERF_BIG_NNZ) ? PERF_BIG_ROUNDS : PERF_SMALL_ROUNDS);

    TVector_SMP in, out;
    DSAFE(vector_create(&in,  matr->cols));
    DSAFE(vector_create(&out, matr->rows));

    for (int i = 0; i < matr->cols; ++i)
        in.val[i] = ((i%2)?(-1):1)*0.01;

    desc_t *desc;

    double dt = timer();
    DSAFE(ker->create(&desc, matr, info));
    dt = timer() - dt;
    p->setup_time = dt;

    double best_dt = 1e100;
    for (int r = 0; r < rounds; ++r) {
        double dt = timer();
        for (int b = 0; b < batch; ++b) {
            DSAFE(ker->compute(desc, &in, &out));
        }
        dt = timer() - dt;
        if (dt < best_dt) best_dt = dt;
    }

    best_dt /= batch;

    p->time = best_dt;
    p->gf   = gflop(matr) / best_dt;

    ker->destroy(desc);

    return ERROR_NO_ERROR;
}
