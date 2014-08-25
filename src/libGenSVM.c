/**
 * @file libGenSVM.c
 * @author Gertjan van den Burg
 * @date August 8, 2013
 * @brief Main functions for the GenSVM algorithm
 *
 * @details
 * The functions in this file are all functions needed
 * to calculate the optimal separation boundaries for 
 * a multiclass classification problem, using the 
 * GenSVM algorithm.
 *
 */

#include <cblas.h>
#include <math.h>

#include "libGenSVM.h"
#include "gensvm.h"
#include "gensvm_matrix.h"

inline double rnd() { return (double) rand()/0x7FFFFFFF; }

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
void gensvm_simplex_gen(long K, double *U)
{
	long i, j;
	for (i=0; i<K; i++) {
		for (j=0; j<K-1; j++) {
			if (i <= j) {
				matrix_set(U, K-1, i, j, -1.0/sqrt(2.0*(j+1)*(j+2)));
			} else if (i == j+1) {
				matrix_set(U, K-1, i, j, sqrt((j+1)/(2.0*(j+2))));
			} else {
				matrix_set(U, K-1, i, j, 0.0);
			}
		}
	}
}

/**
 * @brief Generate the category matrix
 *
 * @details
 * Generate the category matrix R. The category matrix has 1's everywhere
 * except at the column corresponding to the label of instance i, there the
 * element is 0.
 *
 * @param[in,out] 	model 		corresponding GenModel
 * @param[in] 		dataset 	corresponding GenData
 *
 */
void gensvm_category_matrix(struct GenModel *model, struct GenData *dataset)
{
	long i, j;
	long n = model->n;
	long K = model->K;

	for (i=0; i<n; i++) {
		for (j=0; j<K; j++) {
			if (dataset->y[i] != j+1)
				matrix_set(model->R, K, i, j, 1.0);
			else
				matrix_set(model->R, K, i, j, 0.0);
		}
	}
}

/**
 * @brief Generate the simplex difference matrix
 *
 * @details
 * The simplex difference matrix is a 3D matrix which is constructed
 * as follows. For each instance i, the difference vectors between the row of 
 * the simplex matrix corresponding to the class label of instance i and the
 * other rows of the simplex matrix are calculated. These difference vectors 
 * are stored in a matrix, which is one horizontal slice of the 3D matrix.
 *
 * @param[in,out] 	model 	the corresponding GenModel
 * @param[in] 		data 	the corresponding GenData
 *
 */
void gensvm_simplex_diff(struct GenModel *model, struct GenData *data)
{
	long i, j, k;
	double value;

	long n = model->n;
	long K = model->K;

	for (i=0; i<n; i++) {
		for (j=0; j<K-1; j++) {
			for (k=0; k<K; k++) {
				value = matrix_get(model->U, K-1, data->y[i]-1, j);
				value -= matrix_get(model->U, K-1, k, j);
				matrix3_set(model->UU, K-1, K, i, j, k, value);
			}
		}
	}
}

/**
 * @brief Calculate the scalar errors
 * 
 * @details
 * Calculate the scalar errors q based on the current estimate of V, and
 * store these in Q. It is assumed that the memory for Q has already been 
 * allocated. In addition, the matrix ZV is calculated here. It is assigned 
 * to a pre-allocated block of memory, which is passed to this function.
 *
 * @param[in,out] 	model 	the corresponding GenModel
 * @param[in] 		data 	the corresponding GenData
 * @param[in,out] 	ZV 	a pointer to a memory block for ZV. On exit
 * 				this block is updated with the new ZV matrix
 * 				calculated with GenModel::V.
 *
 */
void gensvm_calculate_errors(struct GenModel *model, struct GenData *data,
	       	double *ZV)
{
	long i, j, k;
	double a, value;

	long n = model->n;
	long m = model->m;
	long K = model->K;

	cblas_dgemm(
			CblasRowMajor,
			CblasNoTrans,
			CblasNoTrans,
			n,
			K-1,
			m+1,
			1.0,
			data->Z,
			m+1,
			model->V,
			K-1,
			0.0,
			ZV,
			K-1);

	Memset(model->Q, double, n*K);
	for (i=0; i<n; i++) {
		for (j=0; j<K-1; j++) {
			a = matrix_get(ZV, K-1, i, j);
			for (k=0; k<K; k++) {
				value = a * matrix3_get(model->UU, K-1, K, i,
					       	j, k);
				matrix_add(model->Q, K, i, k, value);
			}
		}
	}
}	

/**
 * @brief Calculate the Huber hinge errors
 *
 * @details
 * For each of the scalar errors in Q the Huber hinge errors are 
 * calculated. The Huber hinge is here defined as 
 * @f[
 * 	h(q) = 
 * 		\begin{dcases}
 * 			1 - q - \frac{\kappa + 1}{2} & \text{if } q \leq -\kappa \\
 * 			\frac{1}{2(\kappa + 1)} ( 1 - q)^2 & \text{if } q \in (-\kappa, 1] \\
 * 			0 & \text{if } q > 1
 * 		\end{dcases}
 * @f]
 *
 * @param[in,out] model 	the corresponding GenModel
 */
void gensvm_calculate_huber(struct GenModel *model) 
{
	long i, j;
	double q, value;

	for (i=0; i<model->n; i++) {
		for (j=0; j<model->K; j++) {
			q = matrix_get(model->Q, model->K, i, j);
			value = 0.0;
			if (q <= -model->kappa) {
				value = 1.0 - q - (model->kappa+1.0)/2.0;
			} else if (q <= 1.0) {
				value = 1.0/(2.0*model->kappa+2.0)*pow(1.0 - q, 2.0);
			}
			matrix_set(model->H, model->K, i, j, value);
		}
	}
}

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
void gensvm_seed_model_V(struct GenModel *from_model,
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
 * @brief Use step doubling
 *
 * @details
 * Step doubling can be used to speed up the maorization algorithm. Instead of 
 * using the value at the minimimum of the majorization function, the value
 * ``opposite'' the majorization point is used. This can essentially cut the
 * number of iterations necessary to reach the minimum in half.
 *
 * @param[in] 	model	GenModel containing the augmented parameters
 */
void gensvm_step_doubling(struct GenModel *model)
{
	long i, j;
	double value;

	long m = model->m;
	long K = model->K;

	for (i=0; i<m+1; i++) {
		for (j=0; j<K-1; j++) {
			matrix_mul(model->V, K-1, i, j, 2.0);
			value = - matrix_get(model->Vbar, K-1, i, j);
			matrix_add(model->V, K-1, i, j, value);
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
	long *groups;
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
			model->rho[i] = ((double) n)/((double) (groups[data->y[i]-1]*K));
	} else {
		fprintf(stderr, "Unknown weight specification.\n");
		exit(1);
	}
}

