#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "unit.h"
#include "rcm.h"

int main(int argc, char** argv) {

	if (argc < 2) {
		printf("usage: %s <matrix-in-file>\n", argv[0]);
		exit(2);
	}

	int flag = 2;
	if (argc == 3) sscanf(argv[2], "%d", &flag);
//	printf("flag = %d\n",flag);

	TMatrix_DCRS matr;
	matrix_load(&matr, argv[1]);
	
	if ( flag & 1 ) matrix_show(&matr, 1);

	rcm(&matr);

	if ( flag & 2 ) matrix_show(&matr, 0);

	return 0;
}

