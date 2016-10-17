
/**
 * @file gensvm_debug.h
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Header file for gensvm_debug.c
 *
 * @details
 * Contains defines useful for debugging.
 *
 */

#ifndef GENSVM_DEBUG_H
#define GENSVM_DEBUG_H

#include "gensvm_print.h"
#include "gensvm_sparse.h"

void gensvm_print_matrix(double *M, long rows, long cols);
void gensvm_print_sparse(struct GenSparse *A);

#endif
