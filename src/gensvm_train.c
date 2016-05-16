/**
 * @file gensvm_train.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Main function for training a GenSVM model.
 *
 */

#include "gensvm_train.h"

/**
 * @brief Utility function for training a GenSVM model
 *
 * @details
 * This function organizes model allocation, kernel preprocessing, instance 
 * weight initialization, and model training. It is the function that should 
 * be used for training a single GenSVM model. Note that optionally a seed 
 * model can be passed to the function to seed the V matrix with. When no such 
 * model is used this parameter should be set to NULL.
 *
 * @param[in] 	model 		a GenModel instance
 * @param[in] 	data 		a GenData instance with the training data
 * @param[in] 	seed_model 	an optional GenModel to seed the V matrix
 *
 */
void gensvm_train(struct GenModel *model, struct GenData *data,
		struct GenModel *seed_model)
{
	// copy dataset parameters to model
	model->n = data->n;
	model->m = data->m;
	model->K = data->K;

	// initialize the V matrix (potentially with a seed model)
	gensvm_init_V(seed_model, model, data);

	// allocate model
	gensvm_allocate_model(model);

	// preprocess kernel
	gensvm_kernel_preprocess(model, data);

	// reallocate model for kernels
	gensvm_reallocate_model(model, data->n, data->r);

	// initialize weights
	gensvm_initialize_weights(data, model);

	// start training
	gensvm_optimize(model, data);
}
