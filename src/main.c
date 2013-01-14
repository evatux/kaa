#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#include "core.h"
#include "rcm.h"

typedef struct {
    const char* matr_in_file;
    const char* portrait_file;
	const char* info_file;
	int		make_original;
	int		make_modified;
    real	threshold;
	real	graph_threshold;
} config_t;

void print_usage_and_exit() {
	printf("usage: rcm <matrix-in-file> [..OPTIONS..]\n\n");

	printf("  %-20s %-10s\t%s\n", "-p|--png", "pattern_name", "specify pattern for graphics ouptut file");
	printf("  %-20s %-10s\t%s\n", "-t|--threshold", "eps", "set what small element is for algorithm");
	printf("  %-20s %-10s\t%s\n", "-g|--graph_threshold", "eps", "set what small element is for png");
	printf("  %-20s %-10s\t%s\n", "-i|--info_file", "file_name", "set output info file");
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
	FILE* inf = (config.info_file == NULL)?stdout:fopen(config.info_file, "w");
	if (inf == NULL) inf = stdout;
	fprintf(inf, "RCM v0.1\nSource file: %s\n\n", config.matr_in_file);

	TMatrix_DCSR matr_mod, matr_orig;
	if (config.make_original) matrix_load(&matr_orig, config.matr_in_file);
	if (config.make_modified) matrix_load(&matr_mod , config.matr_in_file);

	if (config.make_original)	fprintf(inf, "Source matrix: [size: %d], [nonz: %d], [band: %d]\n", matr_orig.size, matr_orig.nonz, matrix_get_band(&matr_orig));
	else						fprintf(inf, "Source matrix: [size: %d], [nonz: %d], [band: %d]\n", matr_mod.size , matr_mod.nonz , matrix_get_band(&matr_mod ));

	// original matrix
	if ( config.portrait_file != NULL ) {
		if (config.make_original) 	matrix_portrait_pattern(&matr_orig, config.portrait_file, "_asrc", config.graph_threshold);
		else						matrix_portrait_pattern(&matr_mod , config.portrait_file, "_asrc", config.graph_threshold);
	}

	if (config.make_original) {
		rcm(&matr_orig, 0);
		fprintf(inf, "Original RCM output: [size: %d], [nonz: %d], [band: %d]\n", matr_orig.size, matr_orig.nonz, matrix_get_band(&matr_orig));
		if ( config.portrait_file != NULL ) matrix_portrait_pattern(&matr_orig, config.portrait_file, "_orcm", config.graph_threshold);
	}
	if (config.make_modified) {
		rcm(&matr_mod , config.threshold);
		fprintf(inf, "Modified RCM output: [size: %d], [nonz: %d], [band: %d]\n", matr_mod.size, matr_mod.nonz, matrix_get_band(&matr_mod));
		if ( config.portrait_file != NULL ) matrix_portrait_pattern(&matr_mod,  config.portrait_file, "_zrcm", config.graph_threshold);
	}

	if (config.make_original) matrix_destroy(&matr_orig);
	if (config.make_modified) matrix_destroy(&matr_mod );
	if (config.info_file != NULL) fclose(inf);

	return 0;
}
