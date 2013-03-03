#include "graphics.h"

#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

// =============================================
//      internal subroutines
// =============================================

/* A coloured pixel. */
typedef struct {
    uint8_t red;
    uint8_t green;
    uint8_t blue;
} pixel_t;

/* A picture. */
typedef struct  {
    pixel_t *pixels;
    size_t width;
    size_t height;
} bitmap_t;

static pixel_t * pixel_at (bitmap_t * bitmap, int x, int y)
{
    return bitmap->pixels + bitmap->width * y + x;
}

/* Scale matrix coordinates to picture coordinates */
static pixel_t * matrix_pixel_at (bitmap_t *bitmap, int row, int col, int matr_size) 
{
    if (matr_size > MAX_PNG_SIZE) { // use scaling
        int x, y;
        x = (int) ( (double)(col * MAX_PNG_SIZE)/(double)matr_size );
        y = (int) ( (double)(row * MAX_PNG_SIZE)/(double)matr_size );
        return bitmap->pixels + bitmap->width * y + x;
    } else {
        return bitmap->pixels + bitmap->width * row + col;
    }
}

/* Write "bitmap" to a PNG file specified by "path"; returns 0 on
   success, non-zero on error. */
static int save_png_to_file (bitmap_t *bitmap, const char *path)
{
    FILE * fp;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    size_t x, y;
    png_byte ** row_pointers = NULL;
    /* "status" contains the return value of this function. At first
       it is set to a value which means 'failure'. When the routine
       has finished its work, it is set to a value which means
       'success'. */
    int status = -1;
    /* The following number is set by trial and error only. I cannot
       see where it it is documented in the libpng manual.
    */
    int pixel_size = 3;
    int depth = 8;

    fp = fopen (path, "wb");
    if (! fp) {
        goto fopen_failed;
    }

    png_ptr = png_create_write_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL) {
        goto png_create_write_struct_failed;
    }

    info_ptr = png_create_info_struct (png_ptr);
    if (info_ptr == NULL) {
        goto png_create_info_struct_failed;
    }

    /* Set up error handling. */
    if (setjmp (png_jmpbuf (png_ptr))) {
        goto png_failure;
    }

    /* Set image attributes. */
    png_set_IHDR (png_ptr,
                  info_ptr,
                  bitmap->width,
                  bitmap->height,
                  depth,
                  PNG_COLOR_TYPE_RGB,
                  PNG_INTERLACE_NONE,
                  PNG_COMPRESSION_TYPE_DEFAULT,
                  PNG_FILTER_TYPE_DEFAULT);

    /* Initialize rows of PNG. */
    row_pointers = png_malloc (png_ptr, bitmap->height * sizeof (png_byte *));
    for (y = 0; y < bitmap->height; ++y) {
        png_byte *row = 
            png_malloc (png_ptr, sizeof (uint8_t) * bitmap->width * pixel_size);
        row_pointers[y] = row;
        for (x = 0; x < bitmap->width; ++x) {
            pixel_t * pixel = pixel_at (bitmap, x, y);
            *row++ = pixel->red;
            *row++ = pixel->green;
            *row++ = pixel->blue;
        }
    }

    /* Write the image data to "fp". */
    png_init_io (png_ptr, fp);
    png_set_rows (png_ptr, info_ptr, row_pointers);
    png_write_png (png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

    /* The routine has successfully written the file, so we set
       "status" to a value which indicates success. */
    status = 0;

    for (y = 0; y < bitmap->height; y++) {
        png_free (png_ptr, row_pointers[y]);
    }
    png_free (png_ptr, row_pointers);

 png_failure:
 png_create_info_struct_failed:
    png_destroy_write_struct (&png_ptr, &info_ptr);
 png_create_write_struct_failed:
    fclose (fp);
 fopen_failed:
    return status;
}

/* Given "value" and "max", the maximum value which we expect "value"
   to take, this returns an integer between 0 and 255 proportional to
   "value" divided by "max". */
static int pix(int value, int max)
{
    if (value < 0)
        return 0;
    return (int) (256.0 *((double) (value)/(double) max));
}

static void put_point(pixel_t *pixel, int color) 
{
    pixel->red   = 255 * (CL_RED   & color);
    pixel->green = 255 * (CL_GREEN & color);
    pixel->blue  = 255 * (CL_BLUE  & color);
}

// =============================================
//      external subroutines
// =============================================

int make_matrix_portrait_color(TMatrix_DCSR *matr, const char *filename, real threshold) 
{
    bitmap_t portrait;
    pixel_t *pixel;
    int x, y;
    int size;

    size = (matr->size > MAX_PNG_SIZE) ? MAX_PNG_SIZE : matr->size;
    portrait.width  = size;
    portrait.height = size;

    portrait.pixels = (pixel_t*) malloc(sizeof(pixel_t) * portrait.width * portrait.height);

    for (y = 0; y < portrait.height; y++) {
        for (x = 0; x < portrait.width; x++) {
            put_point(pixel_at(&portrait, x, y), CL_WHITE);
        }
    }

    int i, j;
    int ci = 0;

    for (i = 0; i < matr->size; i++) {
        for (j = 0; j < matr->size; j++) {
            if ( i == j ) put_point(matrix_pixel_at(&portrait, i, j, matr->size), (FABS(matr->diag[i]) < threshold)?CL_RED:CL_BLACK);
            else {
                for (ci = matr->row_ptr[i]; ci < matr->row_ptr[i+1]; ci++)
                    if ( matr->col_ind[ci] == j ) {
//                      put_point(matrix_pixel_at(&portrait, i, j, matr->size), CL_BLACK);
                        if (FABS(matr->val[ci]) > GRAPH_ZERO_THRESHOLD)
                            put_point(matrix_pixel_at(&portrait, i, j, matr->size), CL_BLACK);
//                      printf("%4.1f ", matr->val[ci]);
                        break;
                    }
            }
        }
    }

    if (save_png_to_file (&portrait, filename) != 0) return ERROR_GRAPHICS;

    return ERROR_NO_ERROR;
}

int make_matrix_portrait(TMatrix_DCSR *matr, const char *filename)
{
    return make_matrix_portrait_color(matr, filename, 0.);
}
