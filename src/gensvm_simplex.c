/**
 * @file gensvm_simplex.c
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Function for generating the simplex matrix
 *
 * @details
 * Contains the function for generating the simplex matrix for a given number
 * of classes.
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

#include "gensvm_simplex.h"

/**
 * @brief Generate matrix of simplex vertex coordinates
 *
 * @details
 * Generate the simplex matrix. Each row of the created matrix contains the 
 * coordinate vector of a single vertex of the K-simplex in K-1 dimensions.  
 * The simplex generated is a special simplex with edges of length 1. The 
 * simplex matrix U of the GenModel must already have been allocated.
 *
 * @param[in,out] 	model 	a GenModel structure
 */
void gensvm_simplex(struct GenModel *model)
{
	long i, j, K = model->K;

	for (i=0; i<K; i++) {
		for (j=0; j<K-1; j++) {
			if (i <= j) {
				matrix_set(model->U, K, K-1, i, j,
					       	-1.0/sqrt(2.0*(j+1)*(j+2)));
			} else if (i == j+1) {
				matrix_set(model->U, K, K-1, i, j,
					       	sqrt((j+1)/(2.0*(j+2))));
			} else {
				matrix_set(model->U, K, K-1, i, j, 0.0);
			}
		}
	}
}

/**
 * @brief Generate the simplex difference matrix
 *
 * @details
 * The simplex difference matrix is a 2D block matrix which is constructed
 * as follows. For each class i, we have a block of K rows and K-1 columns.
 * Each row in the block for class i contains a row vector with the difference
 * of the simplex matrix, U(i, :) - U(j, :).
 *
 * In the paper the notation @f$\boldsymbol{\delta}_{kj}'@f$ is used for the
 * difference vector of @f$\textbf{u}_k' - \textbf{u}_j'@f$, where
 * @f$\textbf{u}_k'@f$ corresponds to row k of @f$\textbf{U}@f$. Due to the
 * indexing in the paper being 1-based and C indexing is 0 based, the vector
 * @f$\boldsymbol{\delta}_{kj}'@f$ corresponds to the row (k-1)*K+(j-1) in the
 * UU matrix generated here.
 *
 * @param[in,out] 	model 	the corresponding GenModel
 *
 */
void gensvm_simplex_diff(struct GenModel *model)
{
	long i, j, l, K = model->K;
	double value;

	// UU is a 2D block matrix, where block i has the differences:
	// U(i, :) - U(j, :) for all j
	for (i=0; i<K; i++) {
		for (j=0; j<K; j++) {
			for (l=0; l<K-1; l++) {
				value = matrix_get(model->U, K, K-1, i, l);
				value -= matrix_get(model->U, K, K-1, j, l);
				matrix_set(model->UU, K*K, K-1, i*K+j, l, 
						value);
			}
		}
	}
}
