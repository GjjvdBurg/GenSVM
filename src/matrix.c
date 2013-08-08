/**
 * @file matrix.c
 * @author Gertjan van den Burg (burg@ese.eur.nl)
 * @date August 8, 2013
 * @brief Functions facilitating matrix access
 *
 * @details
 * The functions contained in this file are used when
 * accessing or writing to matrices. Seperate functions
 * exist of adding and multiplying existing matrix
 * elements, to ensure this is done in place.
 *
 */

#include "matrix.h"

/**
 * @name matrix_set
 * @brief Set element of matrix 
 * @ingroup matrix
 *
 * Row-Major order is used to set a matrix element. Since matrices
 * of type double are most common in MSVMMaj, this function only
 * deals with that type.
 *
 * @param [in] 	M 	matrix to set element of
 * @param [in]  cols 	number of columns of M
 * @param [in] 	i 	row index of element to write to
 * @param [in] 	j 	column index of element to write to
 * @param [out] val 	value to write to specified element of M
 */
void matrix_set(double *M, long cols, long i, long j, double val)
{
	M[i*cols+j] = val;
}

double matrix_get(double *M, long cols, long i, long j)
{
	return M[i*cols+j];
}

void matrix_add(double *M, long cols, long i, long j, double val)
{
	M[i*cols+j] += val;
}

void matrix_mul(double *M, long cols, long i, long j, double val)
{
	M[i*cols+j] *= val;
}

void matrix3_set(double *M, long N2, long N3, long i, long j, 
		long k, double val)
{
	M[k+N3*(j+N2*i)] = val;
}

double matrix3_get(double *M, long N2, long N3, long i, long j,
		long k)
{
	return M[k+N3*(j+N2*i)];
}


void print_matrix(double *M, long rows, long cols)
{
	long i, j;

	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++)
			info("%8.8f ", matrix_get(M, cols, i, j));
		info("\n");
	}
	info("\n");
}

