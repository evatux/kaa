#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#include "core.h"
#include "rcm.h"
#include "cholesky.h"

#define SAFE( f ) \
	do {												\
		int err = f;									\
		if ( ERROR_NO_ERROR != err ) { 					\
			PRINT_ERROR_MESSAGE(err);					\
			fprintf(inf, "error (main): %d\n", err);	\
			if (config.info_file != NULL) fclose(inf);	\
			exit(err);									\
		}												\
	} while(0)

typedef struct {
    const char* matr_in_file;
    const char* portrait_file;
	const char* info_file;
	int		make_original;
	int		make_modified;
    real	threshold;
	real	graph_threshold;
	real	cheps_threshold;
} config_t;

void print_usage_and_exit() {
	printf("usage: rcm <matrix-in-file> [..OPTIONS..]\n\n");

	printf("  %-20s %-10s\t%s\n", "-i|--info_file", "file_name", "set output info file");
	printf("  %-20s %-10s\t%s\n", "-p|--png", "pattern_name", "specify pattern for graphics ouptut file");
	printf("  %-20s %-10s\t%s\n", "-t|--threshold", "eps", "set what small element is for algorithm");
	printf("  %-20s %-10s\t%s\n", "-g|--graph_threshold", "eps", "set what small element is for png");
	printf("  %-20s %-10s\t%s\n", "-c|--cheps_threshold", "eps", "set what small element is for cholesky");
	printf("  %-20s %-10s\t%s\n", "-o|--original_rcm", "", "make only original rcm");
	printf("  %-20s %-10s\t%s\n", "-m|--modified_rcm", "", "make only modified rcm");

	exit(2);
}

void load_config(int argc, char** argv, config_t *config) {
	int cur_opt;

	if (argc < 2) print_usage_and_exit();

	config->matr_in_file	= argv[1];
	config->portrait_file	= NULL;
	config->info_file		= NULL;
	config->threshold		= EPS_THRESHOLD;
	config->graph_threshold = EPS_THRESHOLD;
	config->cheps_threshold	= CHEPS_THRESHOLD;
	config->make_original	= 1;
	config->make_modified	= 1;
	cur_opt = 1;

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
			sscanf(argv[++cur_opt], "%f", &config->threshold);
		} else
		if (!strcmp(argv[cur_opt], "--graph_threshold") || !strcmp(argv[cur_opt], "-g"))
		{
			sscanf(argv[++cur_opt], "%f", &config->graph_threshold);
		} else
		if (!strcmp(argv[cur_opt], "--cheps_threshold") || !strcmp(argv[cur_opt], "-c"))
		{
			sscanf(argv[++cur_opt], "%f", &config->cheps_threshold);
		} else
		if (!strcmp(argv[cur_opt], "--original_rcm") || !strcmp(argv[cur_opt], "-o"))
		{
			config->make_modified = 0;	config->make_original = 1;
		} else
		if (!strcmp(argv[cur_opt], "--modified_rcm") || !strcmp(argv[cur_opt], "-m"))
		{
			config->make_modified = 1;	config->make_original = 0;
		} else
		{
			print_usage_and_exit(argv[0]);
		}
	}
}

void test(config_t* config) {
	TMatrix_Simple simpM;
	TMatrix_DCSR   d1, d2;
	int i;

	matrix_load(&d1, config->matr_in_file);
	matrix_convert_dcsr2simp(&d1, &simpM);
	matrix_convert_simp2dcsr(&simpM, &d2);

#define CMI(X,I) do { if (d1.X != d2.X) { printf("test failed: (%d) %d != %d\n", I, d1.X, d2.X); exit(2); } } while(0)
#define CMF(X,I) do { if (d1.X != d2.X) { printf("test failed: (%d) %f != %f\n", I, d1.X, d2.X); exit(2); } } while(0)

	printf("\nsize: "); CMI(size,0);
	printf("\nnonz: "); CMI(nonz,0);
	printf("\ndiag: "); for (i = 0; i < d1.size; ++i) CMF(diag[i],i);
	printf("\nval:  "); for (i = 0; i < d1.nonz; ++i) CMF( val[i],i);
	printf("\ncol:  "); for (i = 0; i < d1.nonz; ++i) CMI(col_ind[i],i);
	printf("\nrow:  "); for (i = 0; i < d1.size; ++i) CMI(row_ptr[i],i);

	printf("test passed\n");
	exit(0);
}

int main(int argc, char** argv) {
	config_t config;
	load_config(argc, argv, &config);

//	test(&config);

	FILE* inf = (config.info_file == NULL)?stdout:fopen(config.info_file, "w");
	if (inf == NULL) inf = stdout;
	fprintf(inf, "RCM v0.1\nSource file: %s\n\n", config.matr_in_file);

	TMatrix_DCSR matr_src;

	SAFE( matrix_load(&matr_src, config.matr_in_file) ) ;

	fprintf(inf, "Source matrix:       [size: %d], [nonz: %d], [band: %d]\n", matr_src.size, matr_src.nonz, matrix_get_band(&matr_src));
	if ( config.portrait_file != NULL )	matrix_portrait_pattern(&matr_src, config.portrait_file, "_asrc", config.graph_threshold);

	if (config.make_original) {
		fprintf(inf, "\n=============[ Original RCM ]=============\n");

		TMatrix_DCSR A, LD;
		int neps = 0;
		SAFE(	matrix_copy(&matr_src, &A)	);
		SAFE(	rcm(&A, 0)					);

		fprintf(inf, "\tRCM output:      [nonz: %d], [band: %d]\n", A.nonz, matrix_get_band(&A));
		if ( config.portrait_file != NULL ) matrix_portrait_pattern(&A, config.portrait_file, "_orcm", config.graph_threshold);

		SAFE(	cholesky_decomposition(&A, &LD, config.cheps_threshold, &neps) 	);

		fprintf(inf, "\tCholesky output: [nonz: %d], [cheps: %e], [neps: %d]\n", 2*LD.nonz, config.cheps_threshold, neps);
		if ( config.portrait_file != NULL ) matrix_portrait_pattern(&LD, config.portrait_file, "_ochl", config.graph_threshold);

		matrix_destroy(&A);
		matrix_destroy(&LD);
	}
	if (config.make_modified) {
		fprintf(inf, "\n=============[ Modified RCM ]=============\n");

		TMatrix_DCSR A, LD;
		int neps = 0;
		SAFE(	matrix_copy(&matr_src, &A)	);
		SAFE(	rcm(&A, config.threshold)	);

		fprintf(inf, "\tRCM output:      [nonz: %d], [band: %d]\n", A.nonz, matrix_get_band(&A));
		if ( config.portrait_file != NULL ) matrix_portrait_pattern(&A, config.portrait_file, "_zrcm", config.graph_threshold);

		SAFE(	cholesky_decomposition(&A, &LD, config.cheps_threshold, &neps) 	);

		fprintf(inf, "\tCholesky output: [nonz: %d], [cheps: %e], [neps: %d]\n", 2*LD.nonz, config.cheps_threshold, neps);
		if ( config.portrait_file != NULL ) matrix_portrait_pattern(&LD, config.portrait_file, "_zchl", config.graph_threshold);

		matrix_destroy(&A);
		matrix_destroy(&LD);
	}

	matrix_destroy(&matr_src);
	if (config.info_file != NULL) fclose(inf);
	return 0;
}
