#include "core.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include "unit.h"

#ifndef CONFIG_DISABLE_GRAPHICS
    #include "graphics.h"
#endif

int matrix_load(TMatrix_DCSR *matr, const char* filename)
{
	int size, nonz;
	int i;
	real *diag;
	real *val;
	float v;
	int  *col_ind;
	int  *row_ptr;

	FILE* fp = fopen(filename, "r");
	if ( fp == NULL ) return ERROR_FILE_IO; 
	fscanf(fp, "%d %d", &size, &nonz);
#ifdef _DEBUG_LEVEL_1
	printf("[debug(1)]{matrix_load}: size = %d, nonz = %d\n", size, nonz);
#endif
	matr->size = size;
	matr->nonz = nonz;

	diag    = matr->diag    = (real*)malloc(sizeof(real)*size); 
	val     = matr->val     = (real*)malloc(sizeof(real)*nonz);
	col_ind = matr->col_ind = (int* )malloc(sizeof(int )*nonz);
	row_ptr = matr->row_ptr = (int* )malloc(sizeof(int )*(size+1));

	if ( NULL == diag || NULL == val || NULL == col_ind || NULL == row_ptr ) {
		if (   diag) free(   diag);
		if (    val) free(    val);
		if (col_ind) free(col_ind);
		if (row_ptr) free(row_ptr);

		fprintf(stderr, "error [matrix_load]: memory allocation error\n");
		fclose(fp);
		return ERROR_MEMORY_ALLOCATION;
	}

	// LOAD matrix using `Compressed Row Storage` format
	for (i = 0; i < size; i++) fscanf(fp, "%f ", &v), diag[i] = v;	// diagonal matrix items
	for (i = 0; i < nonz; i++) fscanf(fp, "%f ", &v), val[i]  = v;	// nonzero  matrix items
	for (i = 0; i < nonz; i++) fscanf(fp, "%d ", col_ind + i);		// item's column index
	for (i = 0; i < size; i++) fscanf(fp, "%d ", row_ptr + i);		// row ptr offset
	row_ptr[size] = nonz;

	fclose(fp);

	return ERROR_NO_ERROR;
}

int matrix_save(TMatrix_DCSR *matr, const char *filename)
{
	int size, nonz;
    int i;
    FILE* fp = fopen(filename, "w");
    if ( fp == NULL ) return ERROR_FILE_IO;

    size = matr->size;
    nonz = matr->nonz;

    fprintf(fp, "%d %d\n", size, nonz);

	for (i = 0; i < size; i++) fprintf(fp, "%.8e ", matr->diag[i]);       // diagonal matrix items
    fprintf(fp, "\n");
	for (i = 0; i < nonz; i++) fprintf(fp, "%.8e ", matr->val[i]);        // nonzero  matrix items
    fprintf(fp, "\n");
	for (i = 0; i < nonz; i++) fprintf(fp, "%d ", matr->col_ind[i]);    // item's column index
    fprintf(fp, "\n");
	for (i = 0; i < size; i++) fprintf(fp, "%d ", matr->row_ptr[i]);    // row ptr offset
    fprintf(fp, "\n");

	fclose(fp);

	return ERROR_NO_ERROR;
}

int matrix_portrait(TMatrix_DCSR *matr, const char *filename, real threshold)
{
#ifndef CONFIG_DISABLE_GRAPHICS
	return make_matrix_portrait_color(matr, filename, threshold);
#endif
	return ERROR_UNIMPLEMENTED;
}

int matrix_portrait_pattern(TMatrix_DCSR *matr, const char *pattern, const char *suffix, real threshold)
{
	char str[MAX_FILENAME_LENGTH];
	strcpy(str, pattern);
	strcat(str, suffix);
	strcat(str, ".png");
	return matrix_portrait(matr, str, threshold);
}

void matrix_simp_show(TMatrix_Simple *matr)
{
	int i, j;
	int size = matr->size;
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			if ( i == j )
				printf("%4.1f ", matr->val[i*size+j]);
			else
				if ( matr->val[i*size+j] != 0. )
					printf("%4.1f ", matr->val[i*size+j]);
				else
#ifdef PRINT_ZEROES
					printf("%4.1f ", 0.0);
#else
					printf("     ");
#endif
		}
		printf("\n");
	}
	printf("\n");
}

void matrix_show(TMatrix_DCSR *matr, int flag_ordered)
{
	int i, j;
	int ci = 0;
	int flag;

	for (i = 0; i < matr->size; i++) {
		for (j = 0; j < matr->size; j++) {
			if ( i == j )
				printf("%4.1f ", matr->diag[i]);
			else
			if (flag_ordered) {
				if ( matr->col_ind[ci] == j && matr->row_ptr[i] <= ci && ci < matr->row_ptr[i+1] )
					printf("%4.1f ", matr->val[ci++]);
				else
#ifdef PRINT_ZEROES
					printf("%4.1f ", 0.0);
#else
					printf("     ");
#endif
			}
			else {
				flag = 1;
				for (ci = matr->row_ptr[i]; ci < matr->row_ptr[i+1]; ci++)
					if ( matr->col_ind[ci] == j ) {
						flag = 0;
						printf("%4.1f ", matr->val[ci]);
						break;
					}
#ifdef PRINT_ZEROES
					if (flag) printf("%4.1f ", 0.0);
#else
					if (flag) printf("     ");
#endif
			}
		}
		printf("\n");
	}
	printf("\n");
}

void graph_show(TWGraph *gr)
{
	int i;
	printf("===graph===");
	printf("\nwedge:  ");
	for (int i = 0; i < gr->nonz; i++)  printf("%4.2f ", gr->wedge[i]);
	printf("\nadjncy: ");
	for (int i = 0; i < gr->nonz; i++)  printf("%4d ", gr->adjncy[i]);
	printf("\nxadj:   ");
	for (int i = 0; i <= gr->size; i++) printf("%4d ", gr->xadj[i]);
	printf("\nwvert:  ");
	for (int i = 0; i < gr->size; i++)  printf("%4.2f ", gr->wvert[i]);
	printf("\n============\n");
}

void matrix_simp_destroy(TMatrix_Simple *matr)
{
	if (!matr) return;
	if (matr->val) free(matr->val);
}

void matrix_destroy(TMatrix_DCSR *matr)
{
	if (!matr) return;
	if (matr->diag) 	free(matr->diag   );
	if (matr->val ) 	free(matr->val    );
	if (matr->col_ind)	free(matr->col_ind);
	if (matr->row_ptr)	free(matr->row_ptr);
}

void graph_destroy(TWGraph *gr)
{
	if (!gr) return;
	if (gr->wvert )	free(gr->wvert );
	if (gr->wedge )	free(gr->wedge );
	if (gr->adjncy)	free(gr->adjncy);
	if (gr->xadj  )	free(gr->xadj  );
}

int matrix_get_band(TMatrix_DCSR *matr)
{
	int i, ci = 0;
	int band = 0;

	for (i = 0; i < matr->size; ++i) {
		for (ci = matr->row_ptr[i]; ci < matr->row_ptr[i+1]; ++ci) {
			if ( band < (i - matr->col_ind[ci]) ) band = i - matr->col_ind[ci];
		}
	}
	return band;
}

int build_graph(TWGraph *gr, TMatrix_DCSR *matr) 
{
	int i, j, ci;
	real *wvert, *wedge;
	int  *xadj, *adjncy;
	int size;


	wedge  = gr->wedge  = (real*)malloc(sizeof(real)*(matr->nonz));
	adjncy = gr->adjncy = ( int*)malloc(sizeof( int)*(matr->nonz));
	xadj   = gr->xadj   = ( int*)malloc(sizeof( int)*(matr->size+1));
	wvert  = gr->wvert  = (real*)malloc(sizeof(real)*(matr->size));
	size   = gr->size   = matr->size;
			 gr->nonz   = matr->nonz;

	if ( NULL == wedge || NULL == adjncy || NULL == xadj || NULL == wvert) {
		if (  wedge ) free( wedge);
		if ( adjncy ) free(adjncy);
		if (   xadj ) free(  xadj);
		if (  wvert ) free( wvert);

		fprintf(stderr, "error [graph_build]: memory allocation error\n");
		return ERROR_MEMORY_ALLOCATION;
	}	

	for (i = 0; i < size; i++) wvert[i] = matr->diag[i];		// copy diagonal item --> vertices weight

	ci = 0;
	for (i = 0; i < size; i++) {								// run over all rows
		xadj[i] = ci;
		for (j=matr->row_ptr[i]; j<matr->row_ptr[i+1]; j++) {	// run over all non-zero elements in i-th row
//			if (matr->col_ind[j] > i)							// get only right upper triangle part of matrix (U-part)
				adjncy[ci] = matr->col_ind[j];					// make edge
				wedge[ci]  = matr->val[j];						// set  edge weight
				ci++;
		} 
	}
	xadj[i] = ci;												// xadj[size] = nonz

	return ERROR_NO_ERROR;
}

int build_matrix(TWGraph *gr, TMatrix_DCSR *matr, int flag_new) 
{
	int size, nonz;
	int i, j, ci;
	real *diag;
	real *val;
	int  *col_ind;
	int  *row_ptr;

	size = matr->size = gr->size;
	nonz = matr->nonz = gr->xadj[size];

	if ( flag_new ) {
		diag    = matr->diag    = (real*)malloc(sizeof(real)*size); 
		val     = matr->val     = (real*)malloc(sizeof(real)*nonz);
		col_ind = matr->col_ind = (int* )malloc(sizeof(int )*nonz);
		row_ptr = matr->row_ptr = (int* )malloc(sizeof(int )*(size+1));

		if ( NULL == diag || NULL == val || NULL == col_ind || NULL == row_ptr ) {
			if (   diag) free(   diag);
			if (    val) free(    val);
			if (col_ind) free(col_ind);
			if (row_ptr) free(row_ptr);

			fprintf(stderr, "error [build_matrix]: memory allocation error\n");
			return ERROR_MEMORY_ALLOCATION;
		}
	} else {
		diag    = matr->diag;
		val     = matr->val;
		col_ind = matr->col_ind;
		row_ptr = matr->row_ptr;
	}

	for (i = 0; i < size; i++) diag[i] = gr->wvert[i];	// diagonal matrix items

	ci = 0;
	for (i = 0; i < size; i++) {
		row_ptr[i] = ci;
		for (j=gr->xadj[i]; j<gr->xadj[i+1]; j++) {		// run over all non-zero elements in i-th row
				val[ci]     = gr->wedge[j];				// set  edge weight
				col_ind[ci] = gr->adjncy[j];			// make edge
				ci++;
		} 
	}
	row_ptr[size] = nonz;

	return ERROR_NO_ERROR;
}

int matrix_copy(TMatrix_DCSR *src, TMatrix_DCSR *dst)
{
	int size, nonz;
	int i;

	size = dst->size = src->size;
	nonz = dst->nonz = src->nonz;
	dst->diag    = (real*)malloc(sizeof(real)*size);
	dst->val     = (real*)malloc(sizeof(real)*nonz);
	dst->col_ind = ( int*)malloc(sizeof( int)*nonz);
	dst->row_ptr = ( int*)malloc(sizeof( int)*(size+1));

	if ( NULL == dst->diag || NULL == dst->val || NULL == dst->col_ind || NULL == dst->row_ptr ) {
		if (dst->diag   ) free(dst->diag   );
		if (dst->val    ) free(dst->val    );
		if (dst->col_ind) free(dst->col_ind);
		if (dst->row_ptr) free(dst->row_ptr);

		fprintf(stderr, "error [matrix_copy]: memory allocation error\n");
		return ERROR_MEMORY_ALLOCATION;
	}

	for (i = 0; i < size; ++i) dst->diag[i]		= src->diag[i];
	for (i = 0; i < nonz; ++i) dst->val[i]		= src->val[i];
	for (i = 0; i < nonz; ++i) dst->col_ind[i]	= src->col_ind[i];
	for (i = 0; i <= size; ++i) dst->row_ptr[i]	= src->row_ptr[i];

	return ERROR_NO_ERROR;
}

int matrix_convert_simp2dcsr(TMatrix_Simple *src, TMatrix_DCSR *dst)
{
	int i, j, size, nonz;
	size = dst->size = src->size;

	dst->diag = (real*)malloc(sizeof(real)*size);
	if (NULL == dst->diag) {
		fprintf(stderr, "error: [matrix_convert_simp2dcsr]: memory allocation error\n");
		return ERROR_MEMORY_ALLOCATION;
	}
	for (i = 0; i < size; ++i) dst->diag[i] = src->val[i*size + i];

	nonz = 0;
	for (i = 0; i < size; ++i)
		for (j = 0; j < size; ++j)
			if ( (i != j) && (src->val[i*size + j] != 0.) ) nonz++;
	
	dst->nonz    = nonz;
	dst->val     = (real*)malloc(sizeof(real)*nonz);
	dst->col_ind = ( int*)malloc(sizeof( int)*nonz);
	dst->row_ptr = ( int*)malloc(sizeof( int)*(size+1));

	if ( NULL == dst->val || NULL == dst->col_ind || NULL == dst->row_ptr ) {
		if (dst->val    ) free(dst->val    );
		if (dst->col_ind) free(dst->col_ind);
		if (dst->row_ptr) free(dst->row_ptr);

		fprintf(stderr, "error [matrix_convert_simp2dcsr]: memory allocation error\n");
		return ERROR_MEMORY_ALLOCATION;
	}

	nonz = 0;
	dst->row_ptr[0] = 0;
	for (i = 0; i < size; ++i) {
		for (j = 0; j < size; ++j) {
			if ( (i != j) && (src->val[i*size+j] != 0.) ) {
				dst->val[nonz]     = src->val[i*size+j];
				dst->col_ind[nonz] = j;
				nonz++;
			}
		}
		dst->row_ptr[i+1] = nonz;
	}

	return ERROR_NO_ERROR;
}

int matrix_convert_dcsr2simp(TMatrix_DCSR *src, TMatrix_Simple *dst)
{
	int i, j, size;
	size = dst->size = src->size;

	dst->val = (real*)malloc(sizeof(real)*size*size);
	if ( NULL == dst->val ) {
		fprintf(stderr, "error: [matrix_convert_dcsr2simp]: memory allocation error\n");
		return ERROR_MEMORY_ALLOCATION;
	}

	for (i = 0; i < size; ++i) {
		for (j = 0; j < size; ++j) dst->val[i*size+j] = 0.;

		for (j = src->row_ptr[i]; j < src->row_ptr[i+1]; ++j) {
			dst->val[i*size+src->col_ind[j]] = src->val[j];
		}
	}
	for (i = 0; i < size; ++i) dst->val[i*size+i] = src->diag[i];

	return ERROR_NO_ERROR;
}
