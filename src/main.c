#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#include "core.h"
#include "rcm.h"

typedef struct {
    const char* matr_in_file;
    const char* portrait_file;
    real	threshold;
} config_t;

void print_usage_and_exit(const char* prog) {
	printf("usage: %s <matrix-in-file> [-p|--png pattern_name] [-t|--threshold eps] [-o|--original_rcm]\n", prog);
	exit(2);
}

void load_config(int argc, char** argv, config_t *config) {
	int cur_opt;

	if (argc < 2) print_usage_and_exit(argv[0]);

	config->matr_in_file	= argv[1];
	config->portrait_file	= NULL;
	config->threshold	= EPS_THRESHOLD;
	cur_opt = 1;

	while (++cur_opt < argc) {
		if (!strcmp(argv[cur_opt], "--png") || !strcmp(argv[cur_opt], "-p"))
		{
			config->portrait_file = argv[++cur_opt];
		} else
		if (!strcmp(argv[cur_opt], "--threshold") || !strcmp(argv[cur_opt], "-t"))
		{
			sscanf(argv[++cur_opt], "%f", &config->threshold);
		} else
		if (!strcmp(argv[cur_opt], "--original_rcm") || !strcmp(argv[cur_opt], "-o"))
		{
			config->threshold = 0;
		} else
		{
			print_usage_and_exit(argv[0]);
		}
	}
}

int main(int argc, char** argv) {
	config_t config;
	load_config(argc, argv, &config);

	TMatrix_DCSR matr;
	matrix_load(&matr, config.matr_in_file);

	// original matrix
	if ( config.portrait_file != NULL ) matrix_portrait_pattern(&matr, config.portrait_file, "_asrc");

	rcm(&matr, config.threshold);
	if ( config.portrait_file != NULL ) matrix_portrait_pattern(&matr, config.portrait_file, (config.threshold==0)?"_orcm":"_zrcm");

	return 0;
}
