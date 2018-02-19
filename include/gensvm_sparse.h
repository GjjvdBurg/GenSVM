/**
 * @file gensvm_sparse.h
 * @author G.J.J. van den Burg
 * @date 2016-10-11
 * @brief Header file for gensvm_sparse.c
 *
 * @details
 * Contains declarations of the GenSparse structure and related functions.
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

#ifndef GENSVM_SPARSE_H
#define GENSVM_SPARSE_H

// includes
#include "gensvm_globals.h"

// type declarations

/**
 * @brief A structure to represent a sparse matrix in CSR or CSC format
 *
 * @details
 * This structure holds a sparse matrix in the CSR or CSC format. Refer to
 * <a href="https://en.wikipedia.org/wiki/Sparse_matrix">Wikipedia</a> for 
 * more details. The total storage requirement for this structure is:
 *  - CSR: 2*nnz+n_row+1
 *  - CSC: 2*nnz+n_col+1
 * So it only makes sense to use this format if the number of nonzeros is 
 * smaller than
 *  - CSR: @f$(n_{row}(n_{col} - 1) - 1)/2@f$.
 *  - CSC: @f$(n_{col}(n_{row} - 1) - 1)/2@f$.
 *
 * @param nnz 		number of nonzero elements
 * @param n_row 	rows of the matrix
 * @param n_col 	columns of the matrix
 * @param values 	nonzero values (length nnz)
 * @param ix 		cumulative lengths (CSR: n_row+1, CSC: n_col+1)
 * @param jx 		column (CSR) or row (CSC) indices (length nnz)
 */
struct GenSparse {
	enum { CSR, CSC } type;
	///< type of sparse matrix

	long nnz;
	///< number of nonzero elements
	long n_row;
	///< number of rows of the original matrix
	long n_col;
	///< number of columns of the original matrix

	double *values;
	///< array of matrix values, should be of length nnz
	long *ix;
	///< CSR: array of cumulative row lengths, length n_row + 1
	///< CSC: array of cumulative column lengths, length n_col + 1
	long *jx;
	///< CSR: array of colum indices, length nnz
	///< CSC: array of row indices, length nnz
};

struct GenSparse *gensvm_init_sparse(void);

void gensvm_free_sparse(struct GenSparse *sp);
long gensvm_count_nnz(double *A, long rows, long cols);

bool gensvm_nnz_comparison(long nnz, long rows, long cols);
bool gensvm_nnz_comparison_csr(long nnz, long rows, long cols);
bool gensvm_nnz_comparison_csc(long nnz, long rows, long cols);

bool gensvm_could_sparse(double *A, long rows, long cols);
bool gensvm_could_sparse_csr(double *A, long rows, long cols);
bool gensvm_could_sparse_csc(double *A, long rows, long cols);

struct GenSparse *gensvm_dense_to_sparse(double *A, long rows, long cols);
struct GenSparse *gensvm_dense_to_sparse_csr(double *A, long rows, long cols);
struct GenSparse *gensvm_dense_to_sparse_csc(double *A, long rows, long cols);

double *gensvm_sparse_to_dense(struct GenSparse *A);

struct GenSparse *gensvm_sparse_csr_to_csc(struct GenSparse *spA);

#endif
