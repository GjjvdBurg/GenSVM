/**
 * @file gensvm_debug.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Functions facilitating debugging
 *
 * @details
 * Defines functions useful for debugging matrices.
 *
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
