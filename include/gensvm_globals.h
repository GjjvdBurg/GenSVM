/**
 * @file gensvm_globals.h
 * @author G.J.J. van den Burg
 * @date 2014-01-07
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
 * @copyright
 Copyright 2016, G.J.J. van den Burg.

 This file is part of GenSVM.

 GenSVM is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 GenSVM is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with GenSVM. If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef GENSVM_GLOBALS_H
#define GENSVM_GLOBALS_H

#include "gensvm_memory.h"

// all system libraries are included here
#include <cblas.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// ########################### Type definitions ########################### //

/**
 * @brief type of training used in parameter grid search
 */
typedef enum {
	CV=0, /**< cross validation */
	TT=1  /**< data with existing train/test split */
} TrainType;

/**
 * @brief type of kernel used in training
 */
typedef enum {
	K_LINEAR=0, 	/**< Linear kernel */
	K_POLY=1, 	/**< Polynomial kernel */
	K_RBF=2, 	/**< RBF kernel */
	K_SIGMOID=3,  	/**< Sigmoid kernel */
} KernelType;

// ########################### Global constants ########################### //

/**
 * Maximum line length of files that are read into GenSVM.
 */
#ifndef GENSVM_MAX_LINE_LENGTH
  #define GENSVM_MAX_LINE_LENGTH 1024
#endif

// ###################### Min/Max Utility Functions ####################### //

#ifndef MIN_MAX_DEFINE
  /**
   * Flag to check if minimum/maximum macro's are already defined. This can be 
   * useful when linking.
   */
  #define MIN_MAX_DEFINE
  /**
   * Macro for taking the maximum of two arguments.
   */
  #define maximum(a, b) (a) > (b) ? (a) : (b)
  /**
   * Macro for taking the minimum of two arguments.
   */
  #define minimum(a, b) (a) < (b) ? (a) : (b)
#endif

// ####################### Matrix Utility Functions ####################### //

/**
 * Macro for setting a matrix element (RowMajor order)
 */
#define matrix_set(M, cols, i, j, val) M[(i)*(cols)+j] = val

/**
 * Macro for getting a matrix element (RowMajor order)
 */
#define matrix_get(M, cols, i, j) M[(i)*(cols)+j]

/**
 * Macro for adding to a matrix element (RowMajor order)
 */
#define matrix_add(M, cols, i, j, val) M[(i)*(cols)+j] += val

/**
 * Macro for multiplying a matrix element (RowMajor order)
 */
#define matrix_mul(M, cols, i, j, val) M[(i)*(cols)+j] *= val

// ######################### Other Macros ################################# //

// from: http://stackoverflow.com/q/195975/

#define GENSVM_QUOTE(name) #name
#define GENSVM_STRING(macro) GENSVM_QUOTE(macro)

#ifndef VERSION
  #define VERSION 0.0.0
#endif

#define VERSION_STRING GENSVM_STRING(VERSION)


#endif
