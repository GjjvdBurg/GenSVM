/**
 * @file gensvm_sparse.c
 * @author G.J.J. van den Burg
 * @date 2016-10-11
 * @brief Functions for dealing with sparse data matrices
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

#include "gensvm_sparse.h"

/**
 * @brief Initialize a GenSparse structure
 *
 * @details
 * A GenSparse structure is used to hold a sparse data matrix. We work with
 * Compressed Row Storage (CSR) storage, also known as old Yale format.
 *
 * @return 	initialized GenSparse instance
 */
struct GenSparse *gensvm_init_sparse(void)
{
	struct GenSparse *sp = Malloc(struct GenSparse, 1);
	sp->nnz = 0;
	sp->n_row = 0;
	sp->n_col = 0;

	sp->values = NULL;
	sp->ix = NULL;
	sp->jx = NULL;

	return sp;
}

/**
 * @brief Free an allocated GenSparse structure
 *
 * @details
 * Simply free a previously allocated GenSparse structure by freeing all of
 * its components. Finally, the structure itself is freed, and the pointer is
 * set to NULL for safety.
 *
 * @param[in] 	sp 	GenSparse structure to free
 */
void gensvm_free_sparse(struct GenSparse *sp)
{
	free(sp->values);
	free(sp->ix);
	free(sp->jx);
	free(sp);
	sp = NULL;
}

/**
 * @brief Count the number of nonzeros in a matrix
 *
 * @details
 * This is a utility function to count the number of nonzeros in a dense
 * matrix. This is simply done by comparing with 0.0.
 *
 * @param[in] 	A 	a dense matrix (RowMajor order)
 * @param[in] 	rows 	the number of rows of A
 * @param[in] 	cols 	the number of columns of A
 *
 * @return 		the number of nonzeros in A
 */
long gensvm_count_nnz(double *A, long rows, long cols)
{
	long i, nnz = 0;

	for (i=0; i<rows*cols; i++)
		nnz += (A[i] != 0.0) ? 1 : 0;
	return nnz;
}

/**
 * @brief Compare the number of nonzeros is such that sparsity if worth it
 *
 * @details
 * This is a utility function, see gensvm_could_sparse() for more info.
 *
 * @param[in] 	nnz 	number of nonzero elements
 * @param[in] 	rows 	number of rows
 * @param[in] 	cols 	number of columns
 *
 * @return  		whether or not sparsity is worth it
 */
bool gensvm_nnz_comparison(long nnz, long rows, long cols)
{
	#if MAJOR_ORDER == 'r'
	return gensvm_nnz_comparison_csr(nnz, rows, cols);
	#else
	return gensvm_nnz_comparison_csc(nnz, rows, cols);
	#endif
}

bool gensvm_nnz_comparison_csr(long nnz, long rows, long cols)
{
	return (nnz < (rows*(cols-1.0)-1.0)/2.0);
}

bool gensvm_nnz_comparison_csc(long nnz, long rows, long cols)
{
	return (nnz < (cols*(rows-1.0)-1.0)/2.0);
}

/**
 * @brief Check if it is worthwile to convert to a sparse matrix
 *
 * @details
 * It is only worth to convert to a sparse matrix if the amount of sparsity is
 * sufficient. For this to be the case, the number of nonzeros must be
 * smaller than @f$(n_{row}(n_{col} - 1) - 1)/2@f$. This is tested here. If
 * the amount of nonzero entries is small enough, the function returns the
 * number of nonzeros. If it is too big, it returns -1.
 *
 * @sa
 * gensvm_nnz_comparison()
 *
 * @param[in] 	A 	matrix in dense format (RowMajor order)
 * @param[in] 	rows 	number of rows of A
 * @param[in] 	cols 	number of columns of A
 *
 * @return  		whether or not sparsity is worth it
 */

bool gensvm_could_sparse(double *A, long rows, long cols)
{
	#if MAJOR_ORDER == 'r'
	return gensvm_could_sparse_csr(A, rows, cols);
	#else
	return gensvm_could_sparse_csc(A, rows, cols);
	#endif
}

bool gensvm_could_sparse_csr(double *A, long rows, long cols)
{
	long nnz = gensvm_count_nnz(A, rows, cols);
	return gensvm_nnz_comparison_csr(nnz, rows, cols);
}

bool gensvm_could_sparse_csc(double *A, long rows, long cols)
{
	long nnz = gensvm_count_nnz(A, rows, cols);
	return gensvm_nnz_comparison_csc(nnz, rows, cols);
}

/**
 * @brief Convert a dense matrix to a GenSparse structure if advantageous
 *
 * @details
 * This utility function can be used to convert a dense matrix to a sparse
 * matrix in the form of a GenSparse struture. Note that the allocated memory
 * must be freed by the caller. The user should first check whether using a
 * sparse matrix is worth it by calling gensvm_could_sparse().
 *
 * @param[in] 	A 	a dense matrix in RowMajor order
 * @param[in] 	rows 	number of rows of the matrix A
 * @param[in] 	cols 	number of columns of the matrix A
 *
 * @return 		a GenSparse struct
 */
struct GenSparse *gensvm_dense_to_sparse(double *A, long rows, long cols)
{
	#if MAJOR_ORDER == 'r'
	return gensvm_dense_to_sparse_csr(A, rows, cols);
	#else
	return gensvm_dense_to_sparse_csc(A, rows, cols);
	#endif
}

struct GenSparse *gensvm_dense_to_sparse_csr(double *A, long rows, long cols)
{
	double value;
	long row_cnt;
	long i, j, cnt, nnz = 0;
	struct GenSparse *spA = NULL;

	nnz = gensvm_count_nnz(A, rows, cols);

	spA = gensvm_init_sparse();
	spA->type = CSR;

	spA->nnz = nnz;
	spA->n_row = rows;
	spA->n_col = cols;
	spA->values = Calloc(double, nnz);
	spA->ix = Calloc(long, rows+1);
	spA->jx = Calloc(long, nnz);

	cnt = 0;
	spA->ix[0] = 0;
	for (i=0; i<rows; i++) {
		row_cnt = 0;
		for (j=0; j<cols; j++) {
			value = matrix_get(A, rows, cols, i, j);
			if (value != 0) {
				row_cnt++;
				spA->values[cnt] = value;
				spA->jx[cnt] = j;
				cnt++;
			}
		}
		spA->ix[i+1] = spA->ix[i] + row_cnt;
	}

	return spA;
}

struct GenSparse *gensvm_dense_to_sparse_csc(double *A, long rows, long cols)
{
	double value;
	long col_cnt;
	long i, j, cnt, nnz = 0;
	struct GenSparse *spA = NULL;

	nnz = gensvm_count_nnz(A, rows, cols);

	spA = gensvm_init_sparse();

	spA->type = CSC;
	spA->nnz = nnz;
	spA->n_row = rows;
	spA->n_col = cols;
	spA->values = Calloc(double, nnz);
	spA->ix = Calloc(long, cols+1);
	spA->jx = Calloc(long, nnz);

	cnt = 0;
	spA->ix[0] = 0;
	for (j=0; j<cols; j++) {
		col_cnt = 0;
		for (i=0; i<rows; i++) {
			value = matrix_get(A, rows, cols, i, j);
			if (value != 0) {
				col_cnt++;
				spA->values[cnt] = value;
				spA->jx[cnt] = i;
				cnt++;
			}
		}
		spA->ix[j+1] = spA->ix[j] + col_cnt;
	}

	return spA;
}

/**
 * @brief Convert a GenSparse structure to a dense matrix
 *
 * @details
 * This function converts a GenSparse structure back to a normal dense matrix.  
 * Note that the allocated memory must be freed by the caller.
 *
 * @param[in] 	spA 	a GenSparse structure
 *
 * @return 		a dense matrix
 *
 */
double *gensvm_sparse_to_dense(struct GenSparse *spA)
{
	double value;
	long i, j, a, b;
	double *B = Calloc(double, (spA->n_row)*(spA->n_col));

	long end = (spA->type == CSR) ? spA->n_row : spA->n_col;

	for (a=0; a<end; a++) {
		for (b=spA->ix[a]; b<spA->ix[a+1]; b++) {
			if (spA->type == CSR) {
				i = a;
				j = spA->jx[b];
			} else {
				i = spA->jx[b];
				j = a;
			}
			value = spA->values[b];
			matrix_set(B, spA->n_row, spA->n_col, i, j, value);
		}
	}

	return B;
}

/**
 * @brief Convert a CSR sparse matrix to a CSC sparse matrix
 *
 * @param[in] 	spA 	CSR sparse matrix to convert
 *
 * @return  		CSC sparse matrix
 *
 */
struct GenSparse *gensvm_sparse_csr_to_csc(struct GenSparse *spA)
{
	long a, b, i, j, idx;
	double val;

	struct GenSparse *spB = gensvm_init_sparse();

	spB->type = CSC;
	spB->nnz = spA->nnz;
	spB->n_row = spA->n_row;
	spB->n_col = spA->n_col;
	spB->values = Calloc(double, spB->nnz);
	spB->ix = Calloc(double, spB->n_col+1);
	spB->jx = Calloc(double, spB->nnz);

	spB->ix[0] = 0;

	for (a=0; a<spA->n_row; a++) {
		for (b=spA->ix[a]; b<spA->ix[a+1]; b++) {
			i = a;
			j = spA->jx[b];

			spB->ix[j+1]++;
		}
	}
	for (j=0; j<spB->n_col; j++) {
		spB->ix[j+1] += spB->ix[j];
	}

	long *col_cnt = Calloc(long, spB->n_col);

	for (a=0; a<spA->n_row; a++) {
		for (b=spA->ix[a]; b<spA->ix[a+1]; b++) {
			i = a;
			j = spA->jx[b];
			val = spA->values[b];

			idx = spB->ix[j] + col_cnt[j];

			spB->values[idx] = val;
			spB->jx[idx] = i;

			col_cnt[j]++;
		}
	}

	free(col_cnt);

	return spB;
}
