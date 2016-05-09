/**
 * @file gensvm_matrix.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Header file for gensvm_matrix.c
 *
 * @details
 * Contains function declarations for functions useful for dealing with matrices.
 *
 */

#ifndef GENSVM_MATRIX_H
#define GENSVM_MATRIX_H

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

void print_matrix(double *M, long rows, long cols);

#endif
