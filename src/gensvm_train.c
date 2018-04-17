/**
 * @file gensvm_train.c
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Main function for training a GenSVM model.
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
	long real_seed;

	// copy dataset parameters to model
	model->n = data->n;
	model->m = data->m;
	model->K = data->K;

	// allocate model
	gensvm_allocate_model(model);

	// set the random seed
	real_seed = (model->seed == -1) ? time(NULL) : model->seed;
	gensvm_srand(real_seed);

	// preprocess kernel
	gensvm_kernel_preprocess(model, data);

	// reallocate model for kernels
	gensvm_reallocate_model(model, data->n, data->r);

	// initialize the V matrix (potentially with a seed model)
	gensvm_init_V(seed_model, model, data);

	// initialize weights
	gensvm_initialize_weights(data, model);

	// start training
	gensvm_optimize(model, data);
}
