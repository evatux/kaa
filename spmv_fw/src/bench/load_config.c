#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "core.h"
#include "spmv.h"
#include "bench.h"

static void print_usage_and_exit(const char* argv0)
{
    printf("usage: %s <matrix-in-file> [OPTIONS] -- kernel\"info\"\n", argv0);
    printf("usage: %s -ls\tlist known kernels\n\n", argv0);

    printf("OPTIONS:\n");
    #define FORMAT "\t%-8s%-10s\t%s\n"
    printf(FORMAT, "-f",   "",          "florida input matrix format");
    printf(FORMAT, "-Yq",  "",          "check correctness");
    printf(FORMAT, "-b",   "num",       "perf batch");
    printf(FORMAT, "-r",   "num",       "perf rounds");
    printf(FORMAT, "-s",   "num",       "rhs number");
    #undef FORMAT

    exit(2);
}

static void print_kernel_list_and_exit()
{
    for (int i = 0; i < sizeof(kernel_list)/sizeof(kernel_list[0]); ++i)
    {
        spmv_kernel_t ker = kernel_list[i]();
        printf(" %d %s\n", i + 1, ker.name);
    }
    exit(0);
}

static int is_opt(const char *arg, const char *opt1, const char *opt2)
{
    return !strcmp(arg, opt1) ||
        ((opt2 == NULL) ? 0 : !strcmp(arg, opt2));
}

void load_config(int argc, char** argv, config_t *config)
{
    int cur_opt = 1;

    if (argc < 2 || is_opt(argv[cur_opt], "--help", "-h"))
        print_usage_and_exit(argv[0]);

    if (is_opt(argv[cur_opt], "-ls", "--ls"))
        print_kernel_list_and_exit();

    config->matr_in_file    = argv[1];
    config->info_flag       = 0;
    config->fmc_flag        = 0;
    config->correctness     = CORRECTNESS_NO;
    config->batch           = 0;
    config->rounds          = 0;
    config->rhs             = 1;

    for (int i = 0; i < MAX_KER_NUM; ++i) {
        config->kernel[i] = KER_NONE;
        strcpy(config->info[0], "");
    }
    config->ker_num         = 0;

    while (++cur_opt < argc) {
        if (is_opt(argv[cur_opt], "--info", "-i"))
        {
            config->info_flag = 1;
        } else
        if (is_opt(argv[cur_opt], "--florida_file", "-f"))
        {
            config->fmc_flag = 1;
        } else
        if (is_opt(argv[cur_opt], "--correctness", "-Yq"))
        {
            config->correctness = CORRECTNESS_FAST;
        } else
        if (is_opt(argv[cur_opt], "--correctness-user", "-Y"))
        {
            config->correctness = CORRECTNESS_USER;
        } else
        if (is_opt(argv[cur_opt], "--batch", "-b"))
        {
            config->batch = atoi(argv[++cur_opt]);
        } else
        if (is_opt(argv[cur_opt], "--rounds", "-r"))
        {
            config->rounds = atoi(argv[++cur_opt]);
        } else
        if (is_opt(argv[cur_opt], "--rhs", "-s"))
        {
            config->rhs = atoi(argv[++cur_opt]);
        } else
        if (is_opt(argv[cur_opt], "--", NULL))
        {
            ++cur_opt;
            int i;
            for (i = 0; cur_opt < argc; ++cur_opt, ++i) {
                const char *kers = argv[cur_opt];
                config->kernel[i] = kers[0] - '1';
                strcpy(config->info[i], &kers[1]);
            }
            config->ker_num = i;
        } else
        {
            print_usage_and_exit(argv[0]);
        }
    }
}

