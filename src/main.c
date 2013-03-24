#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#include "core.h"
#include "rcm.h"
#include "md.h"
#include "nd.h"
#include "cholesky.h"
#include "solver.h"

#define RCM 1
#define MD  2
#define ND  3

#define SAFE( f ) \
    do {                                                \
        int err = f;                                    \
        if ( ERROR_NO_ERROR != err ) {                  \
            PRINT_ERROR_MESSAGE(err);                   \
            fprintf(inf, "error (main): %d\n", err);    \
            if (config.info_file != NULL) fclose(inf);  \
            exit(err);                                  \
        }                                               \
    } while(0)

typedef struct {
    const char* matr_in_file;
    const char* matr_out_file;
    const char* portrait_file;
    const char* info_file;
    int     make_original;
    int     make_modified;
    int     save_reordering;
    real    threshold;
    real    graph_threshold;
    real    cheps_threshold;
    real    cheps_substitute;
    int     fmc_flag;
    int     algorithm;
    int     separate_png;
} config_t;

void print_usage_and_exit() {
    printf("usage: rcm <matrix-in-file> [matrix-out-pattern] [..OPTIONS..]\n\n");

    printf("  %-20s %-10s\t%s\n", "-i|--info_file", "file_name", "set output info file");
    printf("  %-20s %-10s\t%s\n", "-p|--png", "pattern_name", "specify pattern for graphics ouptut file");
    printf("  %-20s %-10s\t%s\n", "-t|--threshold", "eps", "set what small element is for algorithm");
    printf("  %-20s %-10s\t%s\n", "-g|--graph_threshold", "eps", "set what small element is for png");
    printf("  %-20s %-10s\t%s\n", "-c|--cheps_threshold", "eps", "set what small element is for cholesky");
    printf("  %-20s %-10s\t%s\n", "-h|--cheps_substitute", "eps", "set substitute value for cholesky");
    printf("  %-20s %-10s\t%s\n", "-o|--original_only", "", "make only original algorithm");
    printf("  %-20s %-10s\t%s\n", "-m|--modified_only", "", "make only modified algorithm");
    printf("  %-20s %-10s\t%s\n", "-r|--reordering_save", "", "save reordering matrix in csr format (matrix-out-pattern should be defined)");
    printf("  %-20s %-10s\t%s\n", "-f|--florida_file", "", "matrix-in-file is represented in florida collection format (.mtx)");
    printf("  %-20s %-10s\t%s\n", "-a|--algorithm", "[rcm|md|nd]", "set reordering algorithm (rcm by default)");
    printf("  %-20s %-10s\t%s\n", "-n|--separate_png", "", "make reordering and cholesky pictures separately");

    exit(2);
}

void load_config(int argc, char** argv, config_t *config) {
#ifdef _DEBUG_LEVEL_0
    int i;
    for (i = 0; i < argc; ++i)
        printf("argv[%d] = %s\n", i, argv[i]);
    printf("====================\n");
#endif
    int cur_opt = 1;
    float v;

    if (argc < 2) print_usage_and_exit();

    config->matr_in_file    = argv[1];
    config->matr_out_file   = NULL;
    config->portrait_file   = NULL;
    config->info_file       = NULL;
    config->threshold       = EPS_THRESHOLD;
    config->graph_threshold = EPS_THRESHOLD;
    config->cheps_threshold = CHEPS_THRESHOLD;
    config->cheps_substitute= CHEPS_THRESHOLD;
    config->make_original   = 1;
    config->make_modified   = 1;
    config->save_reordering = 0;
    config->fmc_flag        = 0;
    config->algorithm       = RCM;
    config->separate_png    = 0;

    if ( argc >= 3 && argv[2][0] != '-' )
    {
        config->matr_out_file = argv[2];
        cur_opt++;
    }

    while (++cur_opt < argc) {
        if (!strcmp(argv[cur_opt], "--png") || !strcmp(argv[cur_opt], "-p"))
        {
            config->portrait_file = argv[++cur_opt];
        } else
        if (!strcmp(argv[cur_opt], "--info_file") || !strcmp(argv[cur_opt], "-i"))
        {
            config->info_file = argv[++cur_opt];
        } else
        if (!strcmp(argv[cur_opt], "--threshold") || !strcmp(argv[cur_opt], "-t"))
        {
            sscanf(argv[++cur_opt], "%f", &v);
            config->threshold = v;
        } else
        if (!strcmp(argv[cur_opt], "--graph_threshold") || !strcmp(argv[cur_opt], "-g"))
        {
            sscanf(argv[++cur_opt], "%f", &v);
            config->graph_threshold = v;
        } else
        if (!strcmp(argv[cur_opt], "--cheps_threshold") || !strcmp(argv[cur_opt], "-c"))
        {
            sscanf(argv[++cur_opt], "%f", &v);
            config->cheps_threshold = v;
        } else
        if (!strcmp(argv[cur_opt], "--cheps_substitute") || !strcmp(argv[cur_opt], "-h"))
        {
            sscanf(argv[++cur_opt], "%f", &v);
            config->cheps_substitute = v;
        } else
        if (!strcmp(argv[cur_opt], "--original_only") || !strcmp(argv[cur_opt], "-o"))
        {
            config->make_modified = 0;  config->make_original = 1;
        } else
        if (!strcmp(argv[cur_opt], "--modified_only") || !strcmp(argv[cur_opt], "-m"))
        {
            config->make_modified = 1;  config->make_original = 0;
        } else
        if (!strcmp(argv[cur_opt], "--reordering_save") || !strcmp(argv[cur_opt], "-r"))
        {
            if (!config->matr_out_file) print_usage_and_exit(argv[0]);
            config->save_reordering = 1;
        } else
        if (!strcmp(argv[cur_opt], "--florida_file") || !strcmp(argv[cur_opt], "-f"))
        {
            config->fmc_flag = 1;
        } else
        if (!strcmp(argv[cur_opt], "--separate_png") || !strcmp(argv[cur_opt], "-n"))
        {
            config->separate_png = 1;
        } else
        if (!strcmp(argv[cur_opt], "--algorithm") || !strcmp(argv[cur_opt], "-a"))
        {
            cur_opt++;

            if (!strcmp(argv[cur_opt], "rcm"))
                config->algorithm = RCM;
            else
            if (!strcmp(argv[cur_opt], "md"))
                config->algorithm = MD;
            else
            if (!strcmp(argv[cur_opt], "nd"))
                config->algorithm = ND;
            else
                 print_usage_and_exit(argv[0]);
        } else
        {
            print_usage_and_exit(argv[0]);
        }
    }
}

int main(int argc, char** argv) {
    config_t config;
    load_config(argc, argv, &config);

    char ALG[10];
    char alg[10];
    char oalg[10];
    char zalg[10];
    char output_filename[MAX_FILENAME_LENGTH];
    int (*reorderer) (TMatrix_DCSR*, real);

    if (config.algorithm == RCM)
    {
        strcpy(ALG, "RCM");
        strcpy(alg, "rcm");
        strcpy(oalg, "orcm");
        strcpy(zalg, "zrcm");
        reorderer = rcm;
    } else
    if (config.algorithm == MD)
    {
        strcpy(ALG, "MD");
        strcpy(alg, "md");
        strcpy(oalg, "omd");
        strcpy(zalg, "zmd");
        reorderer = md;
    } else
    if (config.algorithm == ND)
    {
        strcpy(ALG, "ND");
        strcpy(alg, "nd");
        strcpy(oalg, "ond");
        strcpy(zalg, "znd");
        reorderer = nd;
     }

    FILE* inf = (config.info_file == NULL)?stdout:fopen(config.info_file, "w");
    if (inf == NULL) inf = stdout;
    fprintf(inf, "Reordering\nSource file: %s\n", config.matr_in_file);
    fprintf(inf, "Algorithm: %s\n", ALG);
    fprintf(inf, "Threshold: %.2e\nGraph threshold: %.2e\n", config.threshold, config.graph_threshold);
    fprintf(inf, "Cholesky threshold: %.2e\n", config.cheps_threshold);
    fprintf(inf, "Cholesky substitute: %.2e\n\n", config.cheps_substitute);

    TMatrix_DCSR matr_src;

    if (config.fmc_flag == 0)
        SAFE( matrix_load(&matr_src, config.matr_in_file) ) ;
    else
        SAFE( matrix_load_fmc(&matr_src, config.matr_in_file) ) ;


    fprintf(inf, "Source matrix:       [size: %d], [nonz: %d], [band: %d]\n", matr_src.size, matr_src.nonz, matrix_get_band(&matr_src));
    if ( config.portrait_file != NULL ) matrix_portrait_pattern(&matr_src, config.portrait_file, "a", "src", config.graph_threshold);

    if (config.make_original) {
        fprintf(inf, "\n=============[ Original %s ]=============\n", ALG);

        TMatrix_DCSR A, LD, E;
        int neps = 0, pre_neps = 0;
        int *neps_list = NULL;

        SAFE(   matrix_copy(&matr_src, &A)  );
        SAFE(   reorderer(&A, 0)            );

        if (config.save_reordering) {
            sprintf(output_filename, "%s_%s.csr", config.matr_out_file, oalg);
            SAFE(   matrix_save(&A, output_filename)    );
        }

        cholesky_preanalysis(&A, config.cheps_threshold, &pre_neps);
        SAFE(   cholesky_decomposition(&A, &LD, config.cheps_threshold, config.cheps_substitute, &neps, &neps_list)  );

        fprintf(inf, "\tOutput:          [nonz: %d], [band: %d]\n", A.nonz, matrix_get_band(&A));
        fprintf(inf, "\tPre analysis:    [pre_neps: %d]\n", pre_neps);
        fprintf(inf, "\tCholesky output: [nonz: %d], [cheps: %e], [neps: %d]\n", 2*LD.nonz, config.cheps_threshold, neps);
        if ( config.portrait_file != NULL ) {
            if (config.separate_png) {
                matrix_portrait_with_neps_pattern(&A, config.portrait_file, oalg, "",  config.threshold, neps, neps_list);
                matrix_portrait_with_neps_pattern(&LD, config.portrait_file, oalg, "chl", config.cheps_threshold, neps, neps_list);
            } else {
                matrix_portrait_unite_pattern(&A, &LD, config.portrait_file, oalg, config.threshold, neps, neps_list);
            }
        }

        if (config.matr_out_file) {
            sprintf(output_filename, "%s_%s_%s.sim", config.matr_out_file, oalg, "ide");
            SAFE(   make_ident(&A, &LD, &E, output_filename)     );
            if ( config.portrait_file != NULL ) matrix_portrait_with_neps_pattern(&E, config.portrait_file, oalg, "ide", config.cheps_threshold, neps, neps_list);

            fprintf(inf, "\tAlmost Id: [minE: ?], [maxE: ?], [cond: ?]\n");
        }

        matrix_destroy(&A);
        matrix_destroy(&LD);
        if (config.matr_out_file) matrix_destroy(&E);
        if (neps_list) free(neps_list);
    }
    if (config.make_modified) {
        fprintf(inf, "\n=============[ Modified %s ]=============\n", ALG);

        TMatrix_DCSR A, LD, E;
        int neps = 0, pre_neps = 0;
        int* neps_list;

        SAFE( matrix_copy(&matr_src, &A)        );
        SAFE( reorderer(&A, config.threshold)   );

        if (config.save_reordering) {
            sprintf(output_filename, "%s_%s.csr", config.matr_out_file, zalg);
            SAFE(   matrix_save(&A, output_filename)    );
        }

        cholesky_preanalysis(&A, config.cheps_threshold, &pre_neps);
        SAFE(   cholesky_decomposition(&A, &LD, config.cheps_threshold, config.cheps_substitute, &neps, &neps_list)  );

        fprintf(inf, "\tOutput:      [nonz: %d], [band: %d]\n", A.nonz, matrix_get_band(&A));
        fprintf(inf, "\tPre analysis:    [pre_neps: %d]\n", pre_neps);
        fprintf(inf, "\tCholesky output: [nonz: %d], [cheps: %e], [neps: %d]\n", 2*LD.nonz, config.cheps_threshold, neps);
        if ( config.portrait_file != NULL )
        {
            if (config.separate_png) {
                matrix_portrait_with_neps_pattern(&A, config.portrait_file, zalg, "", config.threshold, neps, neps_list);
                matrix_portrait_with_neps_pattern(&LD, config.portrait_file, zalg, "chl", config.cheps_threshold, neps, neps_list);
            } else {
                matrix_portrait_unite_pattern(&A, &LD, config.portrait_file, zalg, config.threshold, neps, neps_list);
            }
        }

        if (config.matr_out_file) {
            sprintf(output_filename, "%s_%s_%s.sim", config.matr_out_file, zalg, "ide");
            SAFE(   make_ident(&A, &LD, &E, output_filename)     );
            if ( config.portrait_file != NULL ) matrix_portrait_with_neps_pattern(&E, config.portrait_file, zalg, "ide", config.cheps_threshold, neps, neps_list);

            fprintf(inf, "\tAlmost Id: [minE: ?], [maxE: ?], [cond: ?]\n");
        }

        matrix_destroy(&A);
        matrix_destroy(&LD);
        if (config.matr_out_file) matrix_destroy(&E);
        if (neps_list) free(neps_list);
    }

    matrix_destroy(&matr_src);
    if (config.info_file != NULL) fclose(inf);
    return 0;
}
