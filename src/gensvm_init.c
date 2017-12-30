/**
 * @file gensvm_init.c
 * @author G.J.J. van den Burg
 * @date 2014-01-07
 * @brief Functions for initializing model and data structures
 * @details
 *
 * This file contains functions for initializing a GenModel instance
 * and a GenData instance. In addition, default values for these
 * structures are defined here (and only here). Functions for allocating
 * memory for the model structure and freeing of the model and data structures
 * are also included.
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

#include "gensvm_init.h"
#include "gensvm_print.h"

/**
 * @brief Seed the matrix V from an existing model or using rand
 *
 * @details
 * The matrix V must be seeded before the main_loop() can start.
 * This can be done by either seeding it with random numbers or
 * using the solution from a previous model on the same dataset
 * as initial seed. The latter option usually allows for a
 * significant improvement in the number of iterations necessary
 * because the seeded model V is closer to the optimal V.
 *
 * When no seed model is supplied, the rows of V are seeded with random 
 * numbers between the inverse of the minimum and the inverse of the maximum 
 * of the corresponding column of Z. This is done to center the product of the 
 * two in the simplex space.
 *
 * @param[in] 		from_model 	GenModel from which to copy V
 * @param[in,out] 	to_model 	GenModel to which V will be copied
 * @param[in] 		data 		GenData structure with the data
 */
void gensvm_init_V(struct GenModel *from_model,
	       	struct GenModel *to_model, struct GenData *data)
{
	long i, j, k, jj_start, jj_end, jj;
	double cmin, cmax, value, rnd;
	double *col_min = NULL,
	       *col_max = NULL;

	long n = data->n;
	long m = data->m;
	long K = data->K;

	if (from_model == NULL) {
		col_min = Calloc(double, m+1);
		col_max = Calloc(double, m+1);
		for (j=0; j<m+1; j++) {
			col_min[j] = 1.0e100;
			col_max[j] = -1.0e100;
		}

		if (data->Z == NULL) {
			// sparse matrix
			long *visit_count = Calloc(long, m+1);
			for (i=0; i<n; i++) {
				jj_start = data->spZ->ia[i];
				jj_end = data->spZ->ia[i+1];
				for (jj=jj_start; jj<jj_end; jj++) {
					j = data->spZ->ja[jj];
					value = data->spZ->values[jj];

					col_min[j] = minimum(col_min[j], value);
					col_max[j] = maximum(col_max[j], value);
					visit_count[j]++;
				}
			}
			// correction in case the minimum or maximum is 0
			for (j=0; j<m+1; j++) {
				if (visit_count[j] < n) {
					col_min[j] = minimum(col_min[j], 0.0);
					col_max[j] = maximum(col_max[j], 0.0);
				}
			}
			free(visit_count);
		} else {
			// dense matrix
			for (i=0; i<n; i++) {
				for (j=0; j<m+1; j++) {
					value = matrix_get(data->Z, m+1, i, j);
					col_min[j] = minimum(col_min[j], value);
					col_max[j] = maximum(col_max[j], value);
				}
			}
		}
		for (j=0; j<m+1; j++) {
			cmin = (fabs(col_min[j]) < 1e-10) ? -1 : col_min[j];
			cmax = (fabs(col_max[j]) < 1e-10) ? 1 : col_max[j];
			for (k=0; k<K-1; k++) {
				rnd = ((double) rand()) / ((double) RAND_MAX);
				value = 1.0/cmin + (1.0/cmax - 1.0/cmin)*rnd;
				matrix_set(to_model->V, K-1, j, k, value);
			}
		}
		free(col_min);
		free(col_max);
	} else {
		for (i=0; i<m+1; i++) {
			for (j=0; j<K-1; j++) {
				value = matrix_get(from_model->V, K-1, i, j);
				matrix_set(to_model->V, K-1, i, j, value);
			}
		}
	}
}

/**
 * @brief Initialize instance weights
 *
 * @details
 * Instance weights can for instance be used to add additional weights to
 * instances of certain classes. Two default weight possibilities are
 * implemented here. The first is unit weights, where each instance gets
 * weight 1.
 *
 * The second are group size correction weights, which are calculated as
 * @f[
 * 	\rho_i = \frac{n}{Kn_k} ,
 * @f]
 * where @f$ n_k @f$ is the number of instances in group @f$ k @f$ and
 * @f$ y_i = k @f$.
 *
 * @param[in] 		data 	GenData with the dataset
 * @param[in,out] 	model 	GenModel with the weight specification. On
 * 				exit GenModel::rho contains the instance
 * 				weights.
 */
void gensvm_initialize_weights(struct GenData *data, struct GenModel *model)
{
	long *groups = NULL;
	long i;

	long n = model->n;
	long K = model->K;

	if (model->weight_idx == 1) {
		for (i=0; i<n; i++)
			model->rho[i] = 1.0;
	}
	else if (model->weight_idx == 2) {
		groups = Calloc(long, K);
		for (i=0; i<n; i++)
			groups[data->y[i]-1]++;
		for (i=0; i<n; i++)
			model->rho[i] = ((double) n)/((double) (
						groups[data->y[i]-1]*K));
	} else {
		// LCOV_EXCL_START
		gensvm_error("[GenSVM Error]: Unknown weight specification.\n");
		exit(EXIT_FAILURE);
		// LCOV_EXCL_STOP
	}
	free(groups);
}
