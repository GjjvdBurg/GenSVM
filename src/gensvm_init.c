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
	model->t = Calloc(double, K-1);
	model->V = Calloc(double, (m+1)*(K-1));
	model->Vbar = Calloc(double, (m+1)*(K-1));
	model->U = Calloc(double, K*(K-1));
	model->UU = Calloc(double, n*K*(K-1));
	model->Q = Calloc(double, n*K);
	model->H = Calloc(double, n*K);
	model->R = Calloc(double, n*K);
	model->rho = Calloc(double, n);
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

	if (model->n == n && model->m == m)
		return;
	if (model->n != n) {
		model->UU = Realloc(model->UU, double, n*K*(K-1));
		Memset(model->UU, double, n*K*(K-1));

		model->Q = Realloc(model->Q, double, n*K);
		Memset(model->Q, double, n*K);

		model->H = Realloc(model->H, double, n*K);
		Memset(model->H, double, n*K);

		model->R = Realloc(model->R, double, n*K);
		Memset(model->R, double, n*K);

		model->rho = Realloc(model->rho, double, n);
		Memset(model->rho, double, n);

		model->n = n;
	}
	if (model->m != m) {
		model->W = Realloc(model->W, double, m*(K-1));
		Memset(model->W, double, m*(K-1));

		model->V = Realloc(model->V, double, (m+1)*(K-1));
		Memset(model->V, double, (m+1)*(K-1));

		model->Vbar = Realloc(model->Vbar, double, (m+1)*(K-1));
		Memset(model->Vbar, double, (m+1)*(K-1));

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
