#include "core.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include "common.h"
#include "unit.h"

#ifndef CONFIG_DISABLE_GRAPHICS
    #include "graphics.h"
#endif

int matrix_create(TMatrix_CSR *matr, int rows, int cols, int nonz, int clean_flag)
{
    matr->rows = rows;
    matr->cols = cols;
    matr->nonz = nonz;
    matr->val     = (real*)malloc(sizeof(real)*nonz);
    matr->col_ind = (int* )malloc(sizeof(int )*nonz);
    matr->row_ptr = (int* )malloc(sizeof(int )*(rows+1));

    if ( NULL == matr->val  ||
         NULL == matr->col_ind || 
         NULL == matr->row_ptr ) 
    {
        if (matr->val)     free(matr->val);
        if (matr->col_ind) free(matr->col_ind);
        if (matr->row_ptr) free(matr->row_ptr);
        matr->rows = matr->cols = matr->nonz = 0;

        D(fprintf(stderr, "error [matrix_create]: memory allocation error\n"));
        return ERROR_MEMORY;
    }

    matr->row_ptr[0] = 0;

    if (clean_flag)
    {
        int i;
        for (i = 0; i < nonz; ++i) matr->val[i]  = 0.;
    }

    return ERROR_NO_ERROR;
}

int matrix_load(TMatrix_CSR *matr, const char* filename)
{
    int rows, cols, nonz;
    int i;
    float v;
    real *val;
    int  *col_ind;
    int  *row_ptr;

    FILE* fp = fopen(filename, "r");
    if ( fp == NULL ) return DE(ERROR_FILE_IO); 
    fscanf(fp, "%d %d %d", &rows, &cols, &nonz);
    DL(1, fprintf(stderr, "[debug] %s: rows = %d, cols = %d, nonz = %d\n",
                filename, rows, cols, nonz));

    matr->rows = rows;
    matr->cols = cols;
    matr->nonz = nonz;

    val     = matr->val     = (real*)malloc(sizeof(real)*nonz);
    col_ind = matr->col_ind = (int* )malloc(sizeof(int )*nonz);
    row_ptr = matr->row_ptr = (int* )malloc(sizeof(int )*(rows+1));

    if ( NULL == val || NULL == col_ind || NULL == row_ptr ) {
        if (    val) free(    val);
        if (col_ind) free(col_ind);
        if (row_ptr) free(row_ptr);

        fclose(fp);
        return DE(ERROR_MEMORY);
    }

    // LOAD matrix using `Compressed Row Storage` format
    for (i = 0; i < nonz; i++) fscanf(fp, "%f ", &v), val[i]  = v;  // nonzero  matrix items
    for (i = 0; i < nonz; i++) fscanf(fp, "%d ", col_ind + i);      // item's column index
    for (i = 0; i < rows; i++) fscanf(fp, "%d ", row_ptr + i);      // row ptr offset
    row_ptr[rows] = nonz;

    fclose(fp);

    return ERROR_NO_ERROR;
}

int matrix_load_fmc(TMatrix_CSR *matr, const char* filename)
{
    int rows, cols, nonz;
    int it, i, j;
    int err;
    float v;
    char str[MAX_FILENAME_LENGTH];
    TMatrix_SMP A;
    FILE* fp = fopen(filename, "r");
    if ( fp == NULL ) return DE(ERROR_FILE_IO);

    while (fgets(str, MAX_FILENAME_LENGTH, fp))
    {
        if (str[0] == '%') continue;
        else break;
    }

    sscanf(str, "%d %d %d", &rows, &cols, &nonz);
    A.rows = rows;
    A.cols = cols;

    A.val = (real*)malloc(sizeof(real)*rows*cols);
    if (A.val == NULL ) {
        fclose(fp);
        return DE(ERROR_MEMORY);
    }

    for (i = 0; i < rows*cols; ++i) A.val[i] = 0.;

    for (it = 0; it < nonz; ++it) {
        fscanf(fp, "%d %d %f", &i, &j, &v);
        i--; j--;
        A.val[i*cols+j] = A.val[j*cols+i] = v;
    }

    fclose(fp);

    err = matrix_convert_smp2csr(&A, matr);
    matrix_smp_destroy(&A);

    return err;
}

int matrix_save(TMatrix_CSR *matr, const char *filename)
{
    int rows, cols, nonz;
    int i;
    FILE* fp = fopen(filename, "w");
    if ( fp == NULL ) return DE(ERROR_FILE_IO);

    rows = matr->rows;
    cols = matr->cols;
    nonz = matr->nonz;

    fprintf(fp, "%d %d %d\n", rows, cols, nonz);

    for (i = 0; i < nonz; i++) fprintf(fp, "%.8e ", matr->val[i]);
    fprintf(fp, "\n");
    for (i = 0; i < nonz; i++) fprintf(fp, "%d ", matr->col_ind[i]);
    fprintf(fp, "\n");
    for (i = 0; i <= rows; i++) fprintf(fp, "%d ", matr->row_ptr[i]);
    fprintf(fp, "\n");

    fclose(fp);

    return ERROR_NO_ERROR;
}

int matrix_portrait(TMatrix_CSR *matr, const char *filename)
{
#ifndef CONFIG_DISABLE_GRAPHICS
    return make_matrix_portrait(matr, filename);
#endif
    return DE(ERROR_UNIMPLEMENTED);
}

void matrix_smp_show(TMatrix_SMP *matr)
{
    int i, j;
    int rows = matr->rows;
    int cols = matr->cols;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            if ( i == j ) {
                printf("%4.1f ", matr->val[i*cols+j]);
            } else
            {
                if ( matr->val[i*cols+j] != 0. )
                {
                    printf("%4.1f ", matr->val[i*cols+j]);
                } else
                {
#ifdef PRINT_ZEROES
                    printf("%4.1f ", 0.0);
#else
                    printf("     ");
#endif
                }
            }
        }
        printf("\n");
    }
    printf("\n");
}

void matrix_show(TMatrix_CSR *matr, int flag_ordered)
{
    int i, j;
    int ci = 0;
    int flag;

    for (i = 0; i < matr->rows; i++) {
        for (j = 0; j < matr->cols; j++) {
            if (flag_ordered) {
                if ( matr->col_ind[ci] == j && matr->row_ptr[i] <= ci && ci < matr->row_ptr[i+1] )
                {
                    printf("%4.1f ", matr->val[ci++]);
                } else
                {
#ifdef PRINT_ZEROES
                    printf("%4.1f ", 0.0);
#else
                    printf("     ");
#endif
                }
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

void matrix_smp_destroy(TMatrix_SMP *matr)
{
    if (!matr) return;
    if (matr->val) free(matr->val);
}

void matrix_destroy(TMatrix_CSR *matr)
{
    if (!matr) return;
    matr->rows = matr->cols = matr->nonz = 0;
    if (matr->val )     free(matr->val    );
    if (matr->col_ind)  free(matr->col_ind);
    if (matr->row_ptr)  free(matr->row_ptr);
}

int matrix_copy(TMatrix_CSR *src, TMatrix_CSR *dst)
{
    int rows, cols, nonz;
    int i;

    rows = dst->rows = src->rows;
    cols = dst->cols = src->cols;
    nonz = dst->nonz = src->nonz;
    dst->val     = (real*)malloc(sizeof(real)*nonz);
    dst->col_ind = ( int*)malloc(sizeof( int)*nonz);
    dst->row_ptr = ( int*)malloc(sizeof( int)*(rows+1));

    if ( NULL == dst->val || NULL == dst->col_ind || NULL == dst->row_ptr ) {
        if (dst->val    ) free(dst->val    );
        if (dst->col_ind) free(dst->col_ind);
        if (dst->row_ptr) free(dst->row_ptr);

        return DE(ERROR_MEMORY);
    }

    for (i = 0; i  < nonz; ++i) dst->val[i]     = src->val[i];
    for (i = 0; i  < nonz; ++i) dst->col_ind[i] = src->col_ind[i];
    for (i = 0; i <= rows; ++i) dst->row_ptr[i] = src->row_ptr[i];

    return ERROR_NO_ERROR;
}

int matrix_convert_smp2csr(TMatrix_SMP *src, TMatrix_CSR *dst)
{
    int i, j, rows, cols, nonz;
    rows = dst->rows = src->rows;
    cols = dst->cols = src->cols;

    nonz = 0;
    for (i = 0; i < rows; ++i)
        for (j = 0; j < cols; ++j)
            if ( (i != j) && (src->val[i*cols + j] != 0.) ) nonz++;

    dst->nonz    = nonz;
    dst->val     = (real*)malloc(sizeof(real)*nonz);
    dst->col_ind = ( int*)malloc(sizeof( int)*nonz);
    dst->row_ptr = ( int*)malloc(sizeof( int)*(rows+1));

    if ( NULL == dst->val || NULL == dst->col_ind || NULL == dst->row_ptr ) {
        if (dst->val    ) free(dst->val    );
        if (dst->col_ind) free(dst->col_ind);
        if (dst->row_ptr) free(dst->row_ptr);

        return DE(ERROR_MEMORY);
    }

    nonz = 0;
    dst->row_ptr[0] = 0;
    for (i = 0; i < rows; ++i) {
        for (j = 0; j < cols; ++j) {
            if ( src->val[i*cols+j] != 0.) {
                dst->val[nonz]     = src->val[i*cols+j];
                dst->col_ind[nonz] = j;
                nonz++;
            }
        }
        dst->row_ptr[i+1] = nonz;
    }

    return ERROR_NO_ERROR;
}

int matrix_convert_csr2smp(TMatrix_CSR *src, TMatrix_SMP *dst)
{
    int i, j, rows, cols;
    rows = dst->rows = src->rows;
    cols = dst->cols = src->cols;

    dst->val = (real*)malloc(sizeof(real)*rows*cols);
    if ( NULL == dst->val ) {
        return DE(ERROR_MEMORY);
    }

    for (i = 0; i < rows; ++i) {
        for (j = 0; j < cols; ++j)
            dst->val[i*cols+j] = 0.;

        for (j = src->row_ptr[i]; j < src->row_ptr[i+1]; ++j) {
            dst->val[i*cols+src->col_ind[j]] = src->val[j];
        }
    }

    return ERROR_NO_ERROR;
}
