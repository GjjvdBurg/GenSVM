/**
 * @file gensvm_debug.c
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Functions facilitating debugging
 *
 * @details
 * Defines functions useful for debugging matrices.
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

#include "gensvm_debug.h"

/**
 * @brief Print a dense matrix
 *
 * @details
 * Debug function to print a matrix
 *
 * @param[in] 	M 	matrix
 * @param[in] 	rows 	number of rows of M
 * @param[in] 	cols 	number of columns of M
 */
void gensvm_print_matrix(double *M, long rows, long cols)
{
	long i, j;

	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			if (j > 0)
				note(" ");
			note("%+6.6f", matrix_get(M, cols, i, j));
		}
		note("\n");
	}
	note("\n");
}

/**
 * @brief Print a sparse matrix
 *
 * @details
 * Debug function to print a GenSparse sparse matrix
 *
 * @param[in] 	A 	a GenSparse matrix to print
 *
 */
void gensvm_print_sparse(struct GenSparse *A)
{
	long i;

	// print matrix dimensions
	note("Sparse Matrix:\n");
	note("\tnnz = %li, rows = %li, cols = %li\n", A->nnz, A->n_row,
			A->n_col);

	// print nonzero values
	note("\tvalues = [ ");
	for (i=0; i<A->nnz; i++) {
		if (i != 0) note(", ");
		note("%f", A->values[i]);
	}
	note(" ]\n");

	// print row indices
	note("\tIA = [ ");
	for (i=0; i<A->n_row+1; i++) {
		if (i != 0) note(", ");
		note("%i", A->ia[i]);
	}
	note(" ]\n");

	// print column indices
	note("\tJA = [ ");
	for (i=0; i<A->nnz; i++) {
		if (i != 0) note(", ");
		note("%i", A->ja[i]);
	}
	note(" ]\n");
}
