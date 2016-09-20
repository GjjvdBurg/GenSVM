/**
 * @file gensvm_init.c
 * @author Gertjan van den Burg
 * @date January 7, 2014
 * @brief Functions for initializing model and data structures
 *
 * @details
 * This file contains functions for initializing a GenModel instance
 * and a GenData instance. In addition, default values for these
 * structures are defined here (and only here). Functions for allocating
 * memory for the model structure and freeing of the model and data structures
 * are also included.
 *
 */

#include "gensvm_init.h"
#include "gensvm_print.h"

inline double rnd() { return (double) rand()/0x7FFFFFFF; }

/**
 * @brief seed the matrix V from an existing model or using rand
 *
 * @details
 * The matrix V must be seeded before the main_loop() can start.
 * This can be done by either seeding it with random numbers or
 * using the solution from a previous model on the same dataset
 * as initial seed. The latter option usually allows for a
 * significant improvement in the number of iterations necessary
 * because the seeded model V is closer to the optimal V.
 *
 * @param[in] 		from_model 	GenModel from which to copy V
 * @param[in,out] 	to_model 	GenModel to which V will be copied
 */
void gensvm_init_V(struct GenModel *from_model,
	       	struct GenModel *to_model, struct GenData *data)
{
	long i, j, k;
	double cmin, cmax, value;

	long n = data->n;
	long m = data->m;
	long K = data->K;

	if (from_model == NULL) {
		for (i=0; i<m+1; i++) {
			cmin = 1e100;
			cmax = -1e100;
			for (k=0; k<n; k++) {
				value = matrix_get(data->Z, m+1, k, i);
				cmin = minimum(cmin, value);
				cmax = maximum(cmax, value);
			}
			for (j=0; j<K-1; j++) {
				cmin = (abs(cmin) < 1e-10) ? -1 : cmin;
				cmax = (abs(cmax) < 1e-10) ? 1 : cmax;
				value = 1.0/cmin + (1.0/cmax - 1.0/cmin)*rnd();
				matrix_set(to_model->V, K-1, i, j, value);
			}
		}
	} else {
		for (i=0; i<m+1; i++)
			for (j=0; j<K-1; j++) {
				value = matrix_get(from_model->V, K-1, i, j);
				matrix_set(to_model->V, K-1, i, j, value);
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
		err("[GenSVM Error]: Unknown weight specification.\n");
		exit(EXIT_FAILURE);
		// LCOV_EXCL_STOP
	}
	free(groups);
}
