/**
 * @file libMSVMMaj.c
 * @author Gertjan van den Burg (burg@ese.eur.nl)
 * @date August 8, 2013
 * @brief Main functions for the MSVMMaj algorithm
 *
 * @details
 * The functions in this file are all functions needed
 * to calculate the optimal separation boundaries for 
 * a multiclass classification problem, using the 
 * MSVMMaj algorithm.
 *
 */

#include <cblas.h>
#include <math.h>

#include "libMSVMMaj.h"
#include "MSVMMaj.h"
#include "matrix.h"

inline double rnd() { return (double) rand()/0x7FFFFFFF; }

/**
 * @name msvmmaj_simplex_gen
 * @brief Generate matrix of simplex vertex coordinates
 * @ingroup libMSVMMaj
 * 
 * Generate the simplex matrix. Each row of the created 
 * matrix contains the coordinate vector of a single 
 * vertex of the K-simplex in K-1 dimensions. The simplex
 * generated is a special simplex with edges of length 1.
 * The simplex matrix U must already have been allocated.
 *
 * @param [in] 		K 	number of classes
 * @param [in,out] 	U 	simplex matrix of size K * (K-1)
 */
void msvmmaj_simplex_gen(long K, double *U)
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

/*!
	Generate the category matrix R. The category matrix has 1's everywhere
	except at the column corresponding to the label of instance i.
*/
void msvmmaj_category_matrix(struct MajModel *model, struct MajData *dataset)
{
	long i, j;
	long n = model->n;
	long K = model->K;

	for (i=0; i<n; i++) {
		for (j=0; j<K; j++) {
			if (dataset->y[i] != j+1) {
				matrix_set(model->R, K, i, j, 1.0);
			}
		}
	}
}

/*!
 * Simplex diff
 */
void msvmmaj_simplex_diff(struct MajModel *model, struct MajData *data)
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

/*!
	Calculate the errors Q based on the current value of V.
	It is assumed that the memory for Q has already been allocated. 
	In addition, the matrix ZV is calculated here. It is assigned to a 
	pre-allocated block of memory, since it would be inefficient to keep
	reassigning this block at every iteration.
*/
void msvmmaj_calculate_errors(struct MajModel *model, struct MajData *data, double *ZV)
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
				value = a * matrix3_get(model->UU, K-1, K, i, j, k);
				matrix_add(model->Q, K, i, k, value);
			}
		}
	}
}	

/*!
	Calculate the Huber hinge errors for each error in the matrix Q.
*/
void msvmmaj_calculate_huber(struct MajModel *model) 
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
 * @name msvmmaj_seed_model_V
 * @brief seed the matrix V from an existing model or using rand
 * @ingroup libMSVMMaj
 *
 * The matrix V must be seeded before the main_loop() can start. 
 * This can be done by either seeding it with random numbers or 
 * using the solution from a previous model on the same dataset
 * as initial seed. The latter option usually allows for a 
 * significant improvement in the number of iterations necessary
 * because the seeded model V is closer to the optimal V.
 *
 * @param [in] 	from_model 	model from which to copy V
 * @param [in,out] to_model 	model to which V will be copied
 */
void msvmmaj_seed_model_V(struct MajModel *from_model, struct MajModel *to_model)
{
	long i, j;
	long m = to_model->m;
	long K = to_model->K;
	double value;

	if (from_model == NULL) {
		for (i=0; i<m+1; i++)
			for (j=0; j<K-1; j++)
				matrix_set(to_model->V, K-1, i, j, -1.0+2.0*rnd());
	} else {
		for (i=0; i<m+1; i++)
			for (j=0; j<K-1; j++) {
				value = matrix_get(from_model->V, K-1, i, j);
				matrix_set(to_model->V, K-1, i, j, value);
			}
	}
}

/*!
 * Step doubling
 */

void msvmmaj_step_doubling(struct MajModel *model)
{
	long i, j;

	long m = model->m;
	long K = model->K;

	for (i=0; i<m+1; i++) {
		for (j=0; j<K-1; j++) {
			matrix_mul(model->V, K-1, i, j, 2.0);
			matrix_add(model->V, K-1, i, j, -matrix_get(model->Vbar, K-1, i, j));
		}
	}
}

/*!
 * initialize_weights
 */

void msvmmaj_initialize_weights(struct MajData *data, struct MajModel *model)
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

