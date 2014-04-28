#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "core.h"
#include "spmv.h"
#include "bench.h"

static int check_correctness(config_t *conf, TMatrix_CSR *matr)
{
    int glob_fail = 0;
    corr_t c;
    c.corr_type = conf->correctness;
    c.matr      = matr;
    c.rhs       = conf->rhs;

    for (int ki = 0; ki < conf->ker_num; ++ki)
    {
        c.err = 0.;
        c.fail = 0;

        spmv_kernel_t ker = kernel_list[conf->kernel[ki]]();
        DSAFE(correctness(&ker, conf->info[ki], &c));

        if (c.fail) glob_fail = 1;
        printf("%-16s%-10.2g%-8s\n", ker.name, c.err,
                (c.fail) ? "FAILED" : "PASSED");
    }
    return !glob_fail;
}

static int check_perf(config_t *conf, TMatrix_CSR *matr)
{
    perf_t p;
    p.matr      = matr;
    p.batch     = conf->batch;
    p.rounds    = conf->rounds;
    p.rhs       = conf->rhs;

    for (int ki = 0; ki < conf->ker_num; ++ki)
    {
        p.setup_time    = 0.;
        p.time          = 0.;
        p.gf            = 0.;

        spmv_kernel_t ker = kernel_list[conf->kernel[ki]]();
        DSAFE(perf(&ker, conf->info[ki], &p));

        printf("%-16s%4d rhs %10.2g s %10.2g s %10.2g gf\n",
                ker.name, p.rhs, p.setup_time, p.time, p.gf);
    }

    return ERROR_NO_ERROR;
}

int main(int argc, char **argv)
{
    config_t conf;
    load_config(argc, argv, &conf);

    TMatrix_CSR _matr;
    TMatrix_CSR *matr = &_matr;
    if (conf.fmc_flag)
        DSAFE(matrix_load_fmc(matr, conf.matr_in_file));
    else
        DSAFE(matrix_load    (matr, conf.matr_in_file));

    if (conf.correctness) return check_correctness(&conf, matr);

    return check_perf(&conf, matr);
}
