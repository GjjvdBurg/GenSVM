/**
 * @file msvmmaj_init.c
 * @author Gertjan van den Burg
 * @date January 7, 2014
 * @brief Functions for initializing model and data structures
 *
 * @details 
 * This file contains functions for initializing a MajModel instance
 * and a MajData instance. In addition, default values for these 
 * structures are defined here (and only here).
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
	model->use_cholesky = false;

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
	data->use_cholesky = false;

	return data;
}	

