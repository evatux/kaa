#ifndef _COMMON_H_
#define _COMMON_H_

/***************************************************
* File: common.h
* This file defines common data types and constants.
* It should be included in every source file
***************************************************/

#include <math.h>

#ifdef  _DEBUG
#define _DEBUG_LEVEL_0
#define _DEBUG_LEVEL_1
#define _DEBUG_LEVEL_2
#define _DEBUG_LEVEL_MD
//#define _DEBUG_LEVEL_CHOLESKY
//#define _DEBUG_LEVEL_SOLVER
#endif

//#define _ZRCM_SMART

#define PRINT_ZEROES
#define MAX_FILENAME_LENGTH 1024

/*
 * Notation: 
 *      Numeration starts from 0 to n-1
 *
 */

// define precision (float by default)
#define _TEST_DOUBLE 
#ifdef _TEST_DOUBLE
typedef double real;
#define FABS fabs
#define SQRT sqrt
#else
typedef float real;
#define FABS fabsf
#define SQRT sqrtf
#endif

#ifdef SGN_WITH_ZERO
#define SGN(x) (((x)<0)?(-1):((x)==0)?(0):(1))
#else
#define SGN(x) (((x)<0)?(-1):(1))
#endif

#define MIN2(x,y) (((x)<(y))?(x):(y))
#define MAX2(x,y) (((x)>(y))?(x):(y))
#define SWAP(_T,__x,__y) do { _T __t = __x; __x = __y; __y = __t; } while(0)

// define ERROR codes
#define ERROR_NO_ERROR          0
#define ERROR_MEMORY_ALLOCATION 2
#define ERROR_FILE_IO           17
#define ERROR_DIV_BY_ZERO       30
#define ERROR_NEGATIVE_SQRT     31
#define ERROR_SMTH_WRONG        50
#define ERROR_GRAPHICS          98
#define ERROR_UNIMPLEMENTED     99

#define PRINT_ERROR_MESSAGE(x) \
    {                                                                                           \
        switch(x) {                                                                             \
            case 2 : fprintf(stderr, "error(%d): memory allocation error\n",    (x)); break;    \
            case 17: fprintf(stderr, "error(%d): file input/output error\n",    (x)); break;    \
            case 30: fprintf(stderr, "error(%d): division by zero\n",           (x)); break;    \
            case 31: fprintf(stderr, "error(%d): sqrt from negative number\n",  (x)); break;    \
            case 50: fprintf(stderr, "error(%d): something goes wrong\n",       (x)); break;    \
            case 98: fprintf(stderr, "error(%d): graphics error\n",             (x)); break;    \
            case 99: fprintf(stderr, "error(%d): unimplemented feature\n",      (x)); break;    \
        }                                                                                       \
    }

#define PRINT_ERROR_MESSAGE_AND_EXIT(x) \
    {                                   \
        PRINT_ERROR_MESSAGE(x);         \
        exit(x);                        \
    }

#endif
