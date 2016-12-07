/**
 * @file gensvm_zv.c
 * @author G.J.J. van den Burg
 * @date 2016-10-17
 * @brief Functions for computing the ZV matrix product
 *
 * @details
 * This file exists because the product Z*V of two matrices occurs both in the
 * computation of the loss function and for predicting class labels. Moreover,
 * a distinction has to be made between dense Z matrices and sparse Z
 * matrices, hence a seperate file is warranted.
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

#include "gensvm_zv.h"

/**
 * @brief Wrapper around sparse/dense versions of this function
 *
 * @details
 * This function tests if the data is stored in dense format or sparse format
 * by testing if GenData::Z is NULL or not, and calls the corresponding
 * version of this function accordingly.
 *
 * @sa
 * gensvm_calculate_ZV_dense(), gensvm_calculate_ZV_sparse()
 *
 * @param[in] 	model 	a GenModel instance holding the model
 * @param[in] 	data 	a GenData instance with the data
 * @param[out]	ZV 	a pre-allocated matrix of appropriate dimensions
 */
void gensvm_calculate_ZV(struct GenModel *model, struct GenData *data,
		double *ZV)
{
	if (data->Z == NULL)
		gensvm_calculate_ZV_sparse(model, data, ZV);
	else
		gensvm_calculate_ZV_dense(model, data, ZV);
}

/**
 * @brief Compute the product Z*V for when Z is a sparse matrix
 *
 * @details
 * This is a simple sparse-dense matrix multiplication, which uses 
 * cblas_daxpy() for each nonzero element of Z, to compute Z*V.
 *
 * @param[in] 	model 	a GenModel instance holding the model
 * @param[in] 	data 	a GenData instance with the data
 * @param[out]	ZV 	a pre-allocated matrix of appropriate dimensions
 */
void gensvm_calculate_ZV_sparse(struct GenModel *model,
		struct GenData *data, double *ZV)
{
	long i, j, jj, jj_start, jj_end, K,
	    n_row = data->spZ->n_row;
	double z_ij;

	K = model->K;

	long *Zia = data->spZ->ia;
	long *Zja = data->spZ->ja;
	double *vals = data->spZ->values;

	for (i=0; i<n_row; i++) {
		jj_start = Zia[i];
		jj_end = Zia[i+1];

		for (jj=jj_start; jj<jj_end; jj++) {
			j = Zja[jj];
			z_ij = vals[jj];

			cblas_daxpy(K-1, z_ij, &model->V[j*(K-1)], 1,
					&ZV[i*(K-1)], 1);
		}
	}
}

/**
 * @brief Compute the product Z*V for when Z is a dense matrix
 *
 * @details
 * This function uses cblas_dgemm() to compute the matrix product between Z 
 * and V.
 *
 * @param[in] 	model 	a GenModel instance holding the model
 * @param[in] 	data 	a GenData instance with the data
 * @param[out]	ZV 	a pre-allocated matrix of appropriate dimensions
 */
void gensvm_calculate_ZV_dense(struct GenModel *model,
		struct GenData *data, double *ZV)
{
	// use n from data, assume m and K are the same between model and data
	long n = data->n;
	long m = model->m;
	long K = model->K;

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, K-1, m+1,
			1.0, data->Z, m+1, model->V, K-1, 0, ZV, K-1);
}
