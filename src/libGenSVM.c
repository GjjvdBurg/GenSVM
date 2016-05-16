/**
 * @file libGenSVM.c
 * @author Gertjan van den Burg
 * @date August 8, 2013
 * @brief Main functions for the GenSVM algorithm
 *
 * @details
 * The functions in this file are all functions needed
 * to calculate the optimal separation boundaries for
 * a multiclass classification problem, using the
 * GenSVM algorithm.
 *
 */

#include <cblas.h>
#include <math.h>

#include "globals.h"
#include "libGenSVM.h"
#include "gensvm.h"
#include "gensvm_matrix.h"


