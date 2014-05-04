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

#define START   1
#define STOP    2

#ifdef __INTEL_COMPILER
long long __rdtsc();

static double timer_rdtsc(int act)
{
    static long long t = 0;
    if (act == START) {
        t = __rdtsc();
        return 0.;
    } else {
        long long t2 = __rdtsc();
        return (t2 - t)/3.3e9;
    }
}
#endif

static double timer_gettimeofday(int act)
{
    static struct timeval tm = { 0 };
    if (act == START)
    {
        gettimeofday( &tm, NULL );
        return 0.;
    } else
    {
        struct timeval tm_fin;
        gettimeofday( &tm_fin, NULL );
        long sec    = tm_fin.tv_sec  - tm.tv_sec;
        long usec   = tm_fin.tv_usec - tm.tv_usec;
        if (usec < 0) {
            sec += 1;
            usec += 1000*1000;
        }
//        printf("%ld %ld\n", sec, usec);
        return (double)sec  + (double)usec / 1e6;
    }
}

static inline double timer(int act)
{
    timer_gettimeofday(act);
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

    timer(START);
    DSAFE(ker->create(&desc, matr, info));
    double dt = timer(STOP);
    p->setup_time = dt;

    double best_dt = 1e100;
    for (int r = 0; r < rounds; ++r) {
        timer(START);
        for (int b = 0; b < batch; ++b) {
            DSAFE(ker->compute(desc, &in, &out));
        }
        dt = timer(STOP);
        if (dt < best_dt) best_dt = dt;
    }

    best_dt /= batch;

    p->time = best_dt;
    p->gf   = gflop(matr) / best_dt;

    ker->destroy(desc);

    return ERROR_NO_ERROR;
}
