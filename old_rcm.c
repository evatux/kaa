#include "rcm.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include "unit.h"

void old_interface_matrix_load(TMatrix_CRS *matr, const char* filename) {
	int size, nonz;
	int i;
	real *val;
	int  *col_ind;
	int  *row_ptr;

	FILE* fp = fopen(filename, "r");
	fscanf(fp, "%d %d", &size, &nonz);
	matr->size = size;
	matr->nonz = nonz;

	val     = matr->val     = (real*)malloc(sizeof(real)*(nonz  ));
	col_ind = matr->col_ind = (int* )malloc(sizeof(int )*(nonz  ));
	row_ptr = matr->row_ptr = (int* )malloc(sizeof(int )*(size+1));

	if ( NULL == val || NULL == col_ind || NULL == row_ptr ) {
		fprintf(stderr, "error [matrix_load]: memory allocation error\n");
		fclose(fp);
		exit(2);
	}

	// LOAD matrix from `Compressed Row Storage` format
	for (i = 0; i < nonz; i++) fscanf(fp, "%f ", val     + i);	// matrix's items
	for (i = 0; i < nonz; i++) fscanf(fp, "%d ", col_ind + i);	// item's column index
	for (i = 0; i < size; i++) fscanf(fp, "%d ", row_ptr + i);	// row ptr offset

	row_ptr[size] = nonz;
	printf("[debug] matr sz[%d] nonz[%d]\n", size, nonz);

	fclose(fp);
}

void old_interface_matrix_show(TMatrix_CRS *matr) {
	int i, j;
	int ci = 0;

	for (i = 0; i < matr->size; i++) {
		for (j = 0; j < matr->size; j++) {
			if ( matr->col_ind[ci] == j && matr->row_ptr[i] <= ci && ci < matr->row_ptr[i+1] )
				printf("%4.1f ", matr->val[ci++]);
			else
				printf("%4.1f ", 0.0);
		}
		printf("\n");
	}
	printf("\n");
}

void old_interface_graph_builder(TGraph *gr, TMatrix_CRS *matr) {
	int i, j;
	int *xadj, *adjncy;
	int x;
	printf("[debug] size+1 = %d, nonz = %d\n", matr->size+1, matr->nonz);
	xadj   = gr->xadj   = (int*)malloc(sizeof(int)*(matr->size+1));
	adjncy = gr->adjncy = (int*)malloc(sizeof(int)*(matr->nonz));
	gr->size = matr->size;

	if ( NULL == xadj || NULL == adjncy ) {
		if ( xadj ) free(xadj);
		fprintf(stderr, "error [graph_build]: memory allocation error\n");
		exit(2);
	}	

	x = 0;
	for (i = 0; i < matr->size; i++) {							// run over all rows
		xadj[i] = x;
		for (j=matr->row_ptr[i]; j<matr->row_ptr[i+1]; j++) {	// run over all non-zero elements in i-th row
//			if (matr->col_ind[j] > i)							// get only right upper triangle part of matrix (U-part)
				adjncy[x++] = matr->col_ind[j];
		} 
	}
	xadj[i] = x;
}

int old_interface_find_periphery(TGraph *gr) {
	int i;
	int level;
	int max_level, max_vertex, min_neighbour;
	int glob_max_level, glob_max_vertex;
	int viewed;
	int xadj_start, xadj_end; 
	int *vertex_level = (int*)malloc(sizeof(int)*(gr->size));

	if ( NULL == vertex_level ) {
		fprintf(stderr, "error [find_periphery]: memory allocation error\n");
		exit(2);
	}

	int g, s;
	int *xadj, *adjncy;
	xadj   = gr->xadj;
	adjncy = gr->adjncy;

	TQueue queue;
	queue_init(&queue);

	max_level  = 1;
	max_vertex = 0;

	do {
		memset(vertex_level, 0, sizeof(int)*(gr->size));
//		for (i = 0; i < gr->size; i++) vertex_level[i] = 0;

		glob_max_level  = max_level;
		glob_max_vertex = max_vertex;

		// Set level for first vertex equals 1
		queue_push(&queue, glob_max_vertex);
		vertex_level[glob_max_vertex] = 1;

		// Set maxmin variables
		max_level     = 1;
		max_vertex    = glob_max_vertex;
		min_neighbour = gr->size;
		viewed = 1;

		while (queue_pop(&queue, &g)) {
			level = vertex_level[g];				// paranet level

			xadj_start = xadj[g];
			xadj_end   = xadj[g+1];

			for (i = xadj_start; i < xadj_end; i++) {
				s = adjncy[i];
				if (!vertex_level[s]) {
					viewed++;						// one more vertex is viewed now
					vertex_level[s] = level + 1;	// increasing level
					queue_push(&queue, s);			// will 
				}
			}

			if (viewed == gr->size) {				// smart thing to decrease number of operations
				if (level > max_level) {			// search for vertex with max level
					max_level     = level;
					max_vertex    = g;
					min_neighbour = xadj_end- xadj_start;
				}
													// but with minimal number of neighbours
				if ( (level == max_level) && (xadj_end - xadj_start < min_neighbour) ) {
					max_vertex = g;
					min_neighbour = xadj_end - xadj_start;
				}
			}
		}
		printf("[debug] current max vertex: %d [level:%d]\n", max_vertex, max_level);
	} while (max_level > glob_max_level);			// repeat until glob_max_level isn't changed

	printf("[debug] glob_max_level: %d\n", glob_max_level);

	free(vertex_level);

	return glob_max_vertex;
}

void old_interface_find_permutation(TGraph *gr, int root, int* perm) {
	// Making permututation vector
	//
	// v-th vertex comes perm[v] 
	// e.g. perm[0] <-- root
	int v, i, s;
	TQueue queue;
	queue_init(&queue);

	int *rperm = (int*)malloc(sizeof(int)*(gr->size));	// but unfortunately we have to use reversed permutation
	if ( NULL == rperm ) {
		fprintf(stderr, "error [find_permutation]: memory allocation error\n");
		exit(2);
	}
	memset(rperm, 0, sizeof(int)*(gr->size));

	// actually it is simple bs
	queue_push(&queue, root);

	s = 0;
	while (queue_pop(&queue, &v)) {
		perm[s++]  = v;
		rperm[v] = 1;									// put connected vertices closer
		for (i = gr->xadj[v]; i < gr->xadj[v+1]; i++)
			if ( !rperm[i] ) queue_push(&queue, i);		// put unindexed vertices to queue 
	}

	free(rperm);
}

//  DEPRECATED
void old_interface_matrix_reorder(TMatrix_CRS *matr, int *perm) {
	int  size         = matr->size;

	real *val_old     = matr->val;
	int  *col_ind_old = matr->col_ind;
	int  *row_ptr_old = matr->row_ptr;

	real *val         = (real*)malloc(sizeof(real)*(matr->nonz  )); 
	int  *col_ind     = ( int*)malloc(sizeof (int)*(matr->nonz  ));
	int  *row_ptr     = ( int*)malloc(sizeof (int)*(matr->size+1));

	if ( NULL == val || NULL == col_ind || NULL == row_ptr ) {
		if ( val     ) free(val);
		if ( col_ind ) free(col_ind);
		fprintf(stderr, "error [matrix_reordering]: memory allocation error\n");
		exit(2);
	}

	int i, i_old;	//  row 
	int j, j_old;	//  offset
	int j0;
	int sw1, sw2;	//  swap indices
	real swap;		//  swap

	j = 0;
	for (i = 0; i < size; i++) {
//		j0 = j;				// offset of the first element in current row
		row_ptr[i] = j;
		i_old = perm[i];

		sw1 = sw2 = j;
		for (j_old = row_ptr_old[i_old]; j_old < row_ptr_old[i_old+1]; j++, j_old++) {
			col_ind[j] = col_ind_old[j_old];
			val[j]     = val_old[j_old];

			if ( col_ind[j] <= i     ) sw1 = j;		// could     be -1
			if ( col_ind[j] == i_old ) sw2 = j;		// shouldn't be -1
		}

		if ( sw1 != sw2 ) {							//			 v-------------v
			if ( col_ind[sw1] == i )				// . . . . [sw1] . . . . [sw2] . . . .
			{
				swap     = val[sw1];
				val[sw1] = val[sw2];
				val[sw2] = swap;
			}										//			 |-------------v
			else if ( col_ind[sw2] < i )	 		// . . . . [sw2] . . . . 00000 . . . .
			{
				swap = val[sw2];
				for ( ; sw2 < sw1; sw2++ ) {
					val[sw2]     = val[sw2+1];
					col_ind[sw2] = col_ind[sw2+1];
				}
				val[sw2]     = swap;
				col_ind[sw2] = i;
			}										//			 v-------------|
			else									// . . . . 00000 . . . . [sw2] . . . .
			{
				swap = val[sw2];
				for ( ; sw2 > sw1; sw2-- ) {
					val[sw2]     = val[sw2-1];
					col_ind[sw2] = col_ind[sw2-1];
				}
				val[sw2]     = swap;
				col_ind[sw2] = i;
			}
//			else
//				fprintf(stderr, "warning [matrix_reoder]: diagonal element is empty");
		}
	}
	row_ptr[size] = matr->nonz;

	free(val_old);
	free(col_ind_old);
	free(row_ptr_old);

	matr->val     = val;
	matr->col_ind = col_ind;
	matr->row_ptr = row_ptr;
}

//	REORDERING: RCM. Reverse Cuthill-McKey
void old_interface_rcm(TMatrix_CRS *matr) {
	TGraph gr;
printf("[debug] graph_builder\n");
	graph_builder(&gr, matr);
printf("[debug] find_periphery\n");
	int root = find_periphery(&gr); 				// RCM starts here (root vertex)

	int *perm = (int*)malloc(sizeof(int)*gr.size); // matrix permutation vector (reodering vector)
	if ( NULL == perm) {
		fprintf(stderr, "error [rcm]: memory allocation error\n");
		exit(2);
	}

printf("[debug] find_permutation\n");
	find_permutation(&gr, root, perm);				// <-- reordering

	free(gr.adjncy);
	free(gr.xadj);

printf("[debug] matrix_reoder\n");
	matrix_reorder(matr, perm);
}

