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

double rnd(void)
{
	return ((double) gensvm_rand()) / 2147483647.0;
}

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
	long a, b, i, j, k, b_start, b_end, end, *visit_count = NULL;
	double cmin, cmax, value;
	double *col_min = NULL,
	       *col_max = NULL;

	// if no model is supplied, or the dimensions of the supplied model 
	// don't match, then we use random initialization.
	if (from_model == NULL || from_model->m != to_model->m || 
			from_model->K != to_model->K) {
		col_min = Calloc(double, to_model->m+1);
		col_max = Calloc(double, to_model->m+1);
		for (j=0; j<to_model->m+1; j++) {
			col_min[j] = 1.0e100;
			col_max[j] = -1.0e100;
		}

		if (data->Z == NULL) {
			// sparse matrix
			visit_count = Calloc(long, to_model->m+1);
			end = ((data->spZ->type == CSR) ? data->spZ->n_row 
					: data->spZ->n_col);
			for (a=0; a<end; a++) {
				b_start = data->spZ->ix[a];
				b_end = data->spZ->ix[a+1];
				for (b=b_start; b<b_end; b++) {
					if (data->spZ->type == CSR) {
						j = data->spZ->jx[b];
					} else {
						j = a;
					}
					value = data->spZ->values[b];
					col_min[j] = minimum(col_min[j], value);
					col_max[j] = maximum(col_max[j], value);
					visit_count[j]++;
				}
			}
			// correction in case the minimum or maximum is 0
			for (j=0; j<to_model->m+1; j++) {
				if (visit_count[j] < data->spZ->n_row) {
					col_min[j] = minimum(col_min[j], 0.0);
					col_max[j] = maximum(col_max[j], 0.0);
				}
			}
			free(visit_count);
		} else {
			// dense matrix
			for (i=0; i<to_model->n; i++) {
				for (j=0; j<to_model->m+1; j++) {
					value = matrix_get(data->Z, to_model->n, 
							to_model->m+1, i, j);
					col_min[j] = minimum(col_min[j], value);
					col_max[j] = maximum(col_max[j], value);
				}
			}
		}
		for (j=0; j<to_model->m+1; j++) {
			cmin = (fabs(col_min[j]) < 1e-10) ? -1 : col_min[j];
			cmax = (fabs(col_max[j]) < 1e-10) ? 1 : col_max[j];
			for (k=0; k<to_model->K-1; k++) {
				value = 1.0/cmin + (1.0/cmax - 1.0/cmin)*rnd();
				matrix_set(to_model->V, to_model->m+1, 
						to_model->K-1, j, k, value);
			}
		}
		free(col_min);
		free(col_max);
	} else {
		for (i=0; i<to_model->m+1; i++) {
			for (j=0; j<to_model->K-1; j++) {
				value = matrix_get(from_model->V, to_model->m+1, 
						to_model->K-1, i, j);
				matrix_set(to_model->V, to_model->m+1, 
						to_model->K-1, i, j, value);
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

	if (model->weight_idx == 0) {
		if (model->rho == NULL) {
			// LCOV_EXCL_START
			gensvm_error("[GenSVM Error]: No raw weights but "
					"weight_idx = 0\n");
			exit(EXIT_FAILURE);
			// LCOV_EXCL_STOP
		}
	}
	else if (model->weight_idx == 1) {
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
		gensvm_error("[GenSVM Error]: Unknown weight specification: "
				"%i.\n", model->weight_idx);
		exit(EXIT_FAILURE);
		// LCOV_EXCL_STOP
	}
	free(groups);
}
