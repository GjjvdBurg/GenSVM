/**
 * @file gensvm_matrix.c
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Functions facilitating matrix access
 *
 * @details
 * The functions contained in this file are used when
 * accessing or writing to matrices. Seperate functions
 * exist of adding and multiplying existing matrix
 * elements, to ensure this is done in place.
 *
 */

#include "gensvm_matrix.h"
#include "gensvm_util.h"

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
void print_matrix(double *M, long rows, long cols)
{
	long i, j;

	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++)
			note("%+6.6f ", matrix_get(M, cols, i, j));
		note("\n");
	}
	note("\n");
}
