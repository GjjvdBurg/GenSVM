/**
 * @file gensvm_simplex.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Function for generating the simplex matrix
 *
 * @details
 * Contains the function for generating the simplex matrix for a given number
 * of classes.
 *
 */

#include "gensvm_simplex.h"

/**
 * @brief Generate matrix of simplex vertex coordinates
 *
 * @details
 * Generate the simplex matrix. Each row of the created
 * matrix contains the coordinate vector of a single
 * vertex of the K-simplex in K-1 dimensions. The simplex
 * generated is a special simplex with edges of length 1.
 * The simplex matrix U must already have been allocated.
 *
 * @param[in] 		K 	number of classes
 * @param[in,out] 	U 	simplex matrix of size K * (K-1)
 */
void gensvm_simplex(long K, double *U)
{
	long i, j;
	for (i=0; i<K; i++) {
		for (j=0; j<K-1; j++) {
			if (i <= j) {
				matrix_set(U, K-1, i, j,
					       	-1.0/sqrt(2.0*(j+1)*(j+2)));
			} else if (i == j+1) {
				matrix_set(U, K-1, i, j,
					       	sqrt((j+1)/(2.0*(j+2))));
			} else {
				matrix_set(U, K-1, i, j, 0.0);
			}
		}
	}
}

