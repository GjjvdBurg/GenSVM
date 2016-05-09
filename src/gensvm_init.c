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

#include <math.h>

#include "gensvm.h"
#include "gensvm_init.h"

/**
 * @brief Initialize a GenModel structure
 *
 * @details
 * A GenModel structure is initialized and the default value for the
 * parameters are set. A pointer to the initialized model is returned.
 *
 * @returns 	initialized GenModel
 */
struct GenModel *gensvm_init_model()
{
	struct GenModel *model = Malloc(struct GenModel, 1);

	// set default values
	model->p = 1.0;
	model->lambda = pow(2, -8.0);
	model->epsilon = 1e-6;
	model->kappa = 0.0;
	model->weight_idx = 1;
	model->kerneltype = K_LINEAR;
	model->kernelparam = NULL;

	model->W = NULL;
	model->t = NULL;
	model->V = NULL;
	model->Vbar = NULL;
	model->U = NULL;
	model->UU = NULL;
	model->Q = NULL;
	model->H = NULL;
	model->R = NULL;
	model->rho = NULL;
	model->data_file = NULL;

	return model;
}

/**
 * @brief Initialize a GenData structure
 *
 * @details
 * A GenData structure is initialized and default values are set.
 * A pointer to the initialized data is returned.
 *
 * @returns 	initialized GenData
 *
 */
struct GenData *gensvm_init_data()
{
	struct GenData *data = Malloc(struct GenData, 1);
	data->Sigma = NULL;
	data->y = NULL;
	data->Z = NULL;
	data->RAW = NULL;

	// set default values
	data->kerneltype = K_LINEAR;
	data->kernelparam = NULL;

	return data;
}

/**
 * @brief Allocate memory for a GenModel
 *
 * @details
 * This function can be used to allocate the memory needed for a GenModel. All
 * arrays in the model are specified and initialized to 0.
 *
 * @param[in] 	model 	GenModel to allocate
 *
 */
void gensvm_allocate_model(struct GenModel *model)
{
	long n = model->n;
	long m = model->m;
	long K = model->K;

	model->W = Calloc(double, m*(K-1));
	if (model->W == NULL) {
		fprintf(stderr, "Failed to allocate memory for W.\n");
		exit(1);
	}

	model->t = Calloc(double, K-1);
	if (model->t == NULL) {
		fprintf(stderr, "Failed to allocate memory for t.\n");
		exit(1);
	}

	model->V = Calloc(double, (m+1)*(K-1));
	if (model->V == NULL) {
		fprintf(stderr, "Failed to allocate memory for V.\n");
		exit(1);
	}

	model->Vbar = Calloc(double, (m+1)*(K-1));
	if (model->Vbar == NULL) {
		fprintf(stderr, "Failed to allocate memory for Vbar.\n");
		exit(1);
	}

	model->U = Calloc(double, K*(K-1));
	if (model->U == NULL) {
		fprintf(stderr, "Failed to allocate memory for U.\n");
		exit(1);
	}

	model->UU = Calloc(double, n*K*(K-1));
	if (model->UU == NULL) {
		fprintf(stderr, "Failed to allocate memory for UU.\n");
		exit(1);
	}

	model->Q = Calloc(double, n*K);
	if (model->Q == NULL) {
		fprintf(stderr, "Failed to allocate memory for Q.\n");
		exit(1);
	}

	model->H = Calloc(double, n*K);
	if (model->H == NULL) {
		fprintf(stderr, "Failed to allocate memory for H.\n");
		exit(1);
	}

	model->R = Calloc(double, n*K);
	if (model->R == NULL) {
		fprintf(stderr, "Failed to allocate memory for R.\n");
		exit(1);
	}

	model->rho = Calloc(double, n);
	if (model->rho == NULL) {
		fprintf(stderr, "Failed to allocate memory for rho.\n");
		exit(1);
	}
}

/**
 * @brief Reallocate memory for GenModel
 *
 * @details
 * This function can be used to reallocate existing memory for a GenModel,
 * upon a change in the model dimensions. This is used in combination with
 * kernels.
 *
 * @param[in] 	model 	GenModel to reallocate
 * @param[in] 	n 	new value of GenModel->n
 * @param[in] 	m 	new value of GenModel->m
 *
 */
void gensvm_reallocate_model(struct GenModel *model, long n, long m)
{
	long K = model->K;
	double *tmp = NULL;

	if (model->n == n && model->m == m)
		return;
	if (model->n != n) {
		tmp = (double *) realloc(model->UU, n*K*(K-1)*sizeof(double));
		if (tmp) {
			Memset(tmp, double, n*K*(K-1));
			model->UU = tmp;
		} else {
			fprintf(stderr, "Failed to reallocate UU\n");
			exit(1);
		}

		tmp = (double *) realloc(model->Q, n*K*sizeof(double));
		if (tmp) {
			Memset(tmp, double, n*K);
			model->Q = tmp;
		} else {
			fprintf(stderr, "Failed to reallocate Q\n");
			exit(1);
		}

		tmp = (double *) realloc(model->H, n*K*sizeof(double));
		if (tmp) {
			Memset(tmp, double, n*K);
			model->H = tmp;
		} else {
			fprintf(stderr, "Failed to reallocate H\n");
			exit(1);
		}

		tmp = (double *) realloc(model->R, n*K*sizeof(double));
		if (tmp) {
			Memset(tmp, double, n*K);
			model->R = tmp;
		} else {
			fprintf(stderr, "Failed to reallocate R\n");
			exit(1);
		}

		tmp = (double *) realloc(model->rho, n*sizeof(double));
		if (tmp) {
			Memset(tmp, double, n);
			model->rho = tmp;
		} else {
			fprintf(stderr, "Failed to reallocte rho\n");
			exit(1);
		}

		model->n = n;
	}
	if (model->m != m) {
		tmp = (double *) realloc(model->W, m*(K-1)*sizeof(double));
		if (tmp) {
			Memset(tmp, double, m*(K-1));
			model->W = tmp;
		} else {
			fprintf(stderr, "Failed to reallocate W\n");
			exit(1);
		}

		tmp = (double *) realloc(model->V, (m+1)*(K-1)*sizeof(double));
		if (tmp) {
			Memset(tmp, double, (m+1)*(K-1));
			model->V = tmp;
		} else {
			fprintf(stderr, "Failed to reallocate V\n");
			exit(1);
		}

		tmp = (double *) realloc(model->Vbar, (m+1)*(K-1)*sizeof(double));
		if (tmp) {
			Memset(tmp, double, (m+1)*(K-1));
			model->Vbar = tmp;
		} else {
			fprintf(stderr, "Failed to reallocate Vbar\n");
			exit(1);
		}

		model->m = m;
	}
}

/**
 * @brief Free allocated GenModel struct
 *
 * @details
 * Simply free a previously allocated GenModel by freeing all its component
 * arrays. Note that the model struct itself is also freed here.
 *
 * @param[in] 	model 	GenModel to free
 *
 */
void gensvm_free_model(struct GenModel *model)
{
	free(model->W);
	free(model->t);
	free(model->V);
	free(model->Vbar);
	free(model->U);
	free(model->UU);
	free(model->Q);
	free(model->H);
	free(model->rho);
	free(model->R);
	free(model->kernelparam);

	free(model);
}

/**
 * @brief Free allocated GenData struct
 *
 * @details
 * Simply free a previously allocated GenData struct by freeing all its
 * components. Note that the data struct itself is also freed here.
 *
 * @param[in] 	data 	GenData struct to free
 *
 */
void gensvm_free_data(struct GenData *data)
{
	if (data == NULL)
		return;

	if (data->Z == data->RAW) {
		free(data->Z);
	} else {
		free(data->Z);
		free(data->RAW);
	}
	free(data->kernelparam);
	free(data->y);
	free(data->Sigma);
	free(data);
}
