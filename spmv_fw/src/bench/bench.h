#ifndef _BENCH_H_
#define _BENCH_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "core.h"
#include "spmv.h"

#define CORRECTNESS_NO      0
#define CORRECTNESS_FAST    1
#define CORRECTNESS_USER    2

#define PERF_BIG_NNZ        (1000*1000)
#define PERF_SMALL_BATCH    20
#define PERF_SMALL_ROUNDS   20
#define PERF_BIG_BATCH      5
#define PERF_BIG_ROUNDS     5

#ifdef  _DOUBLE
#define EPS         1e-12
#define THRESHOLD   EPS
#else
#define EPS         1e-5
#define THRESHOLD   EPS
#endif

static int      correctness_rows        = 100;
static int      correctness_cols        = 200;
static double   correctness_sparsity    = 0.1;

typedef struct {
    const char* matr_in_file;
    int         fmc_flag;
    int         correctness;
    int         batch;          /* default: 0 - auto */
    int         rounds;         /* default: 0 - auto */
    int         rhs;            /* default: 1 */
    int         ker_num;
    int         kernel[MAX_KER_NUM];
    info_t      info[MAX_KER_NUM];
} config_t;

typedef struct {
    int corr_type;
    TMatrix_CSR *matr;
    int rhs;
    double err;         /*out*/
    int fail;           /*out*/
} corr_t;

typedef struct {
    TMatrix_CSR *matr;
    int batch;
    int rounds;
    int rhs;
    double setup_time;  /*out*/
    double time;        /*out*/
    double gf;          /*out*/
} perf_t;

void load_config(int /*argc*/, char** /*argv*/, config_t* /*config*/);
int  correctness(spmv_kernel_t* /*ker*/, info_t /*info*/, corr_t* /*c*/);
int  perf       (spmv_kernel_t* /*ker*/, info_t /*info*/, perf_t* /*p*/);

#endif
