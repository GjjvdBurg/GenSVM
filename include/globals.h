/**
 * @file globals.h
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Global definitions
 *
 * @details
 * This header file contains defines and includes which are used in many
 * parts of the program. Most notably, it includes the gensvm_memory.h header
 * which defines functions for safe memory allocation.
 *
 * Furthermore, a maximum and minimum function are defined here. These
 * functions have their own include guards, to ensure potential linked
 * libraries don't conflict with these definitions.
 *
 */

#ifndef GENSVM_GLOBALS_H
#define GENSVM_GLOBALS_H

#include "gensvm_memory.h"
#include "types.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <cblas.h>

#define MAX_LINE_LENGTH 1024

#ifndef MIN_MAX_DEFINE
#define MIN_MAX_DEFINE
#define maximum(a, b) (a) > (b) ? (a) : (b)
#define minimum(a, b) (a) < (b) ? (a) : (b)
#endif

// Set a matrix element (RowMajor)
#define matrix_set(M, cols, i, j, val) M[(i)*(cols)+j] = val

// Get a matrix element (RowMajor)
#define matrix_get(M, cols, i, j) M[(i)*(cols)+j]

// Add to a matrix element (RowMajor)
#define matrix_add(M, cols, i, j, val) M[(i)*(cols)+j] += val

// Multiply a matrix element (RowMajor)
#define matrix_mul(M, cols, i, j, val) M[(i)*(cols)+j] *= val

// Set a 3D matrix element (N2 = second dim, N3 = third dim, RowMajor)
#define matrix3_set(M, N2, N3, i, j, k, val) M[k+(N3)*(j+(N2)*(i))] = val

// Get a 3D matrix element (N2 = second dim, N3 = third dim, RowMajor)
#define matrix3_get(M, N2, N3, i, j, k) M[k+(N3)*(j+(N2)*(i))]

#endif
