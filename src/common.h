#ifndef _COMMON_H_
#define _COMMON_H_

/***************************************************
* File: common.h
* This file defines common data types and constants.
* It should be included in every source file
***************************************************/

#ifdef  _DEBUG
#define _DEBUG_LEVEL_0
#define _DEBUG_LEVEL_1
#endif

#define PRINT_ZEROES

/*
 * Notation: 
 *		Numeration starts from 0 to n-1
 *
 */

// define precision (float by default)
// #define _TEST_DOUBLE 
#ifdef _TEST_DOUBLE
typedef double real;
#define FABS fabs
#else
typedef float real;
#define FABS fabsf
#endif

// define ERROR codes
#define ERROR_NO_ERROR			0
#define ERROR_MEMORY_ALLOCATION 	2
#define ERROR_FILE_IO 			17

#endif
