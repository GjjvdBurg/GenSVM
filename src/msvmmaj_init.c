/**
 * @file msvmmaj_init.c
 * @author Gertjan van den Burg
 * @date January 7, 2014
 * @brief Functions for initializing model and data structures
 *
 * @details 
 * This file contains functions for initializing a MajModel instance
 * and a MajData instance. In addition, default values for these 
 * structures are defined here (and only here). Functions for allocating
 * memory for the model structure and freeing of the model and data structures
 * are also included.
 *
 */

#include <math.h>

#include "msvmmaj.h"
#include "msvmmaj_init.h"

/**
 * @brief Initialize a MajModel structure
 *
 * @details
 * A MajModel structure is initialized and the default value for the
 * parameters are set. A pointer to the initialized model is returned.
 *
 * @returns 	initialized MajModel
 */
struct MajModel *msvmmaj_init_model()
{
	struct MajModel *model = Malloc(struct MajModel, 1);
	
	// set default values
	model->p = 1.0;
	model->lambda = pow(2, -8.0);
	model->epsilon = 1e-6;
	model->kappa = 0.0;
	model->weight_idx = 1;
	model->kerneltype = K_LINEAR;
	model->kernelparam = NULL;

	return model;
}

/**
 * @brief Initialize a MajData structure
 * 
 * @details
 * A MajData structure is initialized and default values are set. 
 * A pointer to the initialized data is returned.
 *
 * @returns 	initialized MajData
 *
 */
struct MajData *msvmmaj_init_data()
{
	struct MajData *data = Malloc(struct MajData, 1);

	// set default values
	data->kerneltype = K_LINEAR;
	data->kernelparam = NULL;

	return data;
}	

/**
 * @brief Allocate memory for a MajModel
 *
 * @details
 * This function can be used to allocate the memory needed for a MajModel. All
 * arrays in the model are specified and initialized to 0.
 *
 * @param[in] 	model 	MajModel to allocate
 *
 */
void msvmmaj_allocate_model(struct MajModel *model)
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
 * @brief Free allocated MajModel struct
 *
 * @details
 * Simply free a previously allocated MajModel by freeing all its component
 * arrays. Note that the model struct itself is also freed here.
 *
 * @param[in] 	model 	MajModel to free
 *
 */
void msvmmaj_free_model(struct MajModel *model)
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

	free(model);
}

/**
 * @brief Free allocated MajData struct 
 *
 * @details
 * Simply free a previously allocated MajData struct by freeing all its
 * components. Note that the data struct itself is also freed here.
 *
 * @param[in] 	data 	MajData struct to free
 *
 */
void msvmmaj_free_data(struct MajData *data)
{
	free(data->Z);
	free(data->y);
	free(data);
}
