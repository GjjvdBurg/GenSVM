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
 * @brief print a matrix
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
