/**
 * @file gensvm_base.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Functions for initializing GenModel and GenData structures
 *
 * @details
 * This file contains functions for initializing, freeing, allocating, and 
 * reallocating a GenModel instance. It also contains functions for 
 * initializing and freeing a GenData structure. In addition, default values 
 * for these structures are defined here (and only here).
 *
 */

#include "gensvm_base.h"

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
	data = NULL;
}

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

	model->V = NULL;
	model->Vbar = NULL;
	model->U = NULL;
	model->UU = NULL;
	model->Q = NULL;
	model->H = NULL;
	model->rho = NULL;
	model->data_file = NULL;

	return model;
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

	model->V = Calloc(double, (m+1)*(K-1));
	model->Vbar = Calloc(double, (m+1)*(K-1));
	model->U = Calloc(double, K*(K-1));
	model->UU = Calloc(double, K*K*(K-1));
	model->Q = Calloc(double, n*K);
	model->H = Calloc(double, n*K);
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
		model->Q = Realloc(model->Q, double, n*K);
		Memset(model->Q, double, n*K);

		model->H = Realloc(model->H, double, n*K);
		Memset(model->H, double, n*K);

		model->rho = Realloc(model->rho, double, n);
		Memset(model->rho, double, n);

		model->n = n;
	}
	if (model->m != m) {
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
	if (model == NULL)
		return;

	free(model->V);
	free(model->Vbar);
	free(model->U);
	free(model->UU);
	free(model->Q);
	free(model->H);
	free(model->rho);
	free(model->kernelparam);
	free(model->data_file);

	free(model);
	model = NULL;
}

/**
 * @brief Initialize the workspace structure
 *
 * @details
 * During the computations in gensvm_get_update(), a number of matrices need 
 * to be allocated which are used to store intermediate results. To avoid 
 * reallocating these matrices at every call to gensvm_get_update(), we create 
 * here a workspace with these matrices. This workspace is created once by 
 * gensvm_optimize(), and is passed to gensvm_get_update() and 
 * gensvm_get_loss(). See the documentation of the GenWork structure for 
 * information on each allocated field.
 *
 * @param[in] 	model 	a GenModel with the dimensionality of the problem
 * @returns 		an allocated GenWork instance
 *
 */
struct GenWork *gensvm_init_work(struct GenModel *model)
{
	long n = model->n;
	long m = model->m;
	long K = model->K;

	struct GenWork *work = Malloc(struct GenWork, 1);
	work->n = n;
	work->m = m;
	work->K = K;

	work->LZ = Calloc(double, n*(m+1));
	work->ZB = Calloc(double, (m+1)*(K-1)),
	work->ZBc = Calloc(double, (m+1)*(K-1)),
	work->ZAZ = Calloc(double, (m+1)*(m+1)),
	work->ZV = Calloc(double, n*(K-1));
	work->beta = Calloc(double, K-1);

	return work;
}

/**
 * @brief Free an allocated GenWork instance
 *
 * @details
 * This function simply frees every matrix allocated for in the GenWork 
 * workspace.
 *
 * @param[in] 	work 	a pointer to an allocated GenWork instance
 *
 */
void gensvm_free_work(struct GenWork *work)
{
	free(work->LZ);
	free(work->ZB);
	free(work->ZBc);
	free(work->ZAZ);
	free(work->ZV);
	free(work->beta);
	free(work);
	work = NULL;
}

/**
 * @brief Reset all matrices of a GenWork instance
 *
 * @details
 * In the gensvm_get_update() function, it is expected for some matrices that 
 * all their entries are initialized to 0. This function sets the memory for 
 * each of the matrices in the GenWork workspace to 0 to facilitate this.
 *
 * @param[in,out] 	work 	an initialized GenWork instance. On exit, the
 * 				matrices of the GenWork instance are all 
 * 				initialized to 0.
 */
void gensvm_reset_work(struct GenWork *work)
{
	long n = work->n;
	long m = work->m;
	long K = work->K;

	Memset(work->LZ, double, n*(m+1));
	Memset(work->ZB, double, (m+1)*(K-1)),
	Memset(work->ZBc, double, (m+1)*(K-1)),
	Memset(work->ZAZ, double, (m+1)*(m+1)),
	Memset(work->ZV, double, n*(K-1));
	Memset(work->beta, double, K-1);
}
