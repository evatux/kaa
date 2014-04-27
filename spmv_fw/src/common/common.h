#ifndef _COMMON_H_
#define _COMMON_H_

/***************************************************
* File: common.h
* This file defines common data types and constants.
* It should be included in every source file
***************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PRINT_ZEROES
#define MAX_FILENAME_LENGTH 1024

//#define _DOUBLE
#ifdef  _DOUBLE
typedef double real;
#define FABS fabs
#define SQRT sqrt
#else
typedef float real;
#define FABS fabsf
#define SQRT sqrtf
#endif

#ifdef  SGN_WITH_ZERO
#define SGNerr ((err<0)?(-1):(err==0)?(0):(1))
#else
#define SGNerr ((err<0)?(-1):(1))
#endif

#define MIN2(x,y) (((x)<(y))?(x):(y))
#define MAX2(x,y) (((x)>(y))?(x):(y))
#define SWAP(x,y) do { typeof(x) tmp = (x); (x) = (y); (y) = tmp; } while(0)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// define ERROR codes
#define ERROR_NO_ERROR           0
#define ERROR_MEMORY             2
#define ERROR_INVALID_CONF       3
#define ERROR_INVALID_DESC       4
#define ERROR_FILE_IO           17
#define ERROR_DIV_BY_ZERO       30
#define ERROR_NEGATIVE_SQRT     31
#define ERROR_SMTH_WRONG        50
#define ERROR_GRAPHICS          98
#define ERROR_UNIMPLEMENTED     99

static int print_error_message(int err, const char* func)
{
    if (err == ERROR_NO_ERROR) return err;

    fprintf(stderr, "error(%s:%d): ", func, err);
    switch(err) {
        case ERROR_MEMORY:
            fprintf(stderr, "memory allocation error\n");
            break;
        case ERROR_INVALID_CONF:
            fprintf(stderr, "invalid configuration\n");
            break;
        case ERROR_INVALID_DESC:
            fprintf(stderr, "invalid descriptor\n");
            break;
        case ERROR_FILE_IO:
            fprintf(stderr, "file input/output error\n");
            break;
        case ERROR_DIV_BY_ZERO:
            fprintf(stderr, "division by zero\n");
            break;
        case ERROR_NEGATIVE_SQRT:
            fprintf(stderr, "sqrt from negative number\n");
            break;
        case ERROR_SMTH_WRONG:
            fprintf(stderr, "something goes wrong\n");
            break;
        case ERROR_GRAPHICS:
            fprintf(stderr, "graphics error\n");
            break;
        case ERROR_UNIMPLEMENTED:
            fprintf(stderr, "unimplemented feature\n");
            break;
        default:
            fprintf(stderr, "uknown error\n");
    }
    return err;
}

//#define _DEBUG
#define _DEBUG_LEVEL 10
#ifdef  _DEBUG
#define DE(err) print_error_message(err, __func__)
#define DSAFE(x) do { \
    int _err = x; \
    if (_err) { fprintf(stderr, "%s -> err:%d\n", #x, _err); return DE(_err); }\
    } while(0)
#define DL(level,x) do { \
    if (_DEBUG_LEVEL >= level) { \
    fprintf(stderr, "%s(%d:%s)\n", __FILE__, __LINE__, __func__); \
    x; }} while(0)
#define D(x) DL(0,x)
#else
#define DE(err) err
#define DSAFE(x) do { \
    int _err = x; \
    if (_err) return DE(_err);\
    } while(0)
#define DL(level,x)
#define D(x)
#endif

#endif
