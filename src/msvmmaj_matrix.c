/**
 * @file msvmmaj_matrix.c
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

#include "msvmmaj_matrix.h"
#include "util.h"

/**
 * @brief Set element of matrix 
 *
 * @details
 * Row-Major order is used to set a matrix element. Since matrices
 * of type double are most common in MSVMMaj, this function only
 * deals with that type.
 *
 * @param[in] 	M 	matrix to set element of
 * @param[in] 	cols 	number of columns of M
 * @param[in] 	i 	row index of element to write to
 * @param[in] 	j 	column index of element to write to
 * @param[in] 	val 	value to write to specified element of M
 */
void matrix_set(double *M, long cols, long i, long j, double val)
{
	M[i*cols+j] = val;
}

/**
 * @brief Retrieve value from matrix
 * 
 * @details
 * Return a value from a matrix using row-major order.
 *
 * @param[in] 	M 	matrix to retrieve value from
 * @param[in] 	cols 	number of columns of M
 * @param[in] 	i 	row index (starting from 0)
 * @param[in] 	j 	column index (starting from 0)
 * @return 		matrix element at (i, j)
 */
double matrix_get(double *M, long cols, long i, long j)
{
	return M[i*cols+j];
}

/**
 * @brief Add value to matrix element
 *
 * @details
 * This function is added to efficiently add values to matrix
 * elements, without having to use get and set methods.
 *
 * @param[in] 	M 	matrix
 * @param[in] 	cols 	number of columns of M
 * @param[in] 	i  	row index (starting from 0)
 * @param[in] 	j 	column index (starting from 0)
 * @param[in] 	val 	value to add to matrix element (i, j)
 */
void matrix_add(double *M, long cols, long i, long j, double val)
{
	M[i*cols+j] += val;
}

/**
 * @brief Multiply matrix element by value
 * 
 * @details
 * This function is added to efficiently multiply a matrix element
 * by a certain value, without having to use get and set methods.
 *
 * @param[in] 	M 	matrix
 * @param[in] 	cols 	number of columns of M
 * @param[in] 	i 	row index (starting from 0)
 * @param[in] 	j 	column index (starting from 0)
 * @param[in] 	val 	value to multiply matrix element (i, j) with
 */
void matrix_mul(double *M, long cols, long i, long j, double val)
{
	M[i*cols+j] *= val;
}

/**
 * @brief Set element of 3D matrix
 *
 * @details
 * Set an element of a 3D matrix using row-major order. 
 *
 * @param[in] 	M 	matrix
 * @param[in] 	N2 	second dimension of M
 * @param[in] 	N3 	third dimension of M
 * @param[in] 	i 	index along first dimension
 * @param[in] 	j 	index along second dimension
 * @param[in] 	k 	index along third dimension
 * @param[in] 	val 	value to set element (i, j, k) to
 *
 * See:
 * http://en.wikipedia.org/wiki/Row-major_order
 */
void matrix3_set(double *M, long N2, long N3, long i, long j, 
		long k, double val)
{
	M[k+N3*(j+N2*i)] = val;
}

/**
 * @brief Get element of 3D matrix
 *
 * @details
 * Retrieve an element from a 3D matrix.
 *
 * @param[in] 	M 	matrix
 * @param[in] 	N2 	second dimension of M
 * @param[in] 	N3 	third dimension of M
 * @param[in] 	i 	index along first dimension
 * @param[in] 	j 	index along second dimension
 * @param[in] 	k 	index along third dimension
 * @returns 		value at the (i, j, k) element of M
 */
double matrix3_get(double *M, long N2, long N3, long i, long j,
		long k)
{
	return M[k+N3*(j+N2*i)];
}

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
			note("%8.8f ", matrix_get(M, cols, i, j));
		note("\n");
	}
	note("\n");
}
