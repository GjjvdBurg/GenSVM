/**
 * @file gensvm_task.c
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Functions for initializing and freeing a GenTask
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

#include "gensvm_task.h"

/**
 * @brief Initialize a GenTask structure
 *
 * @details
 * A GenTask structure is initialized and the default value for the
 * parameters are set. A pointer to the initialized GenTask is returned.
 *
 * @returns 	initialized GenTask
 */
struct GenTask *gensvm_init_task(void)
{
	struct GenTask *t = Malloc(struct GenTask, 1);

	// for training
	t->ID = -1;
	t->folds = 10;
	t->train_data = NULL;
	t->test_data = NULL;

	// for the model
	t->kerneltype = K_LINEAR;
	t->weight_idx = 1;
	t->p = 1.0;
	t->kappa = 0.0;
	t->lambda = 1.0;
	t->epsilon = 1e-6;
	t->gamma = 1.0;
	t->coef = 0.0;
	t->degree = 2.0;
	t->max_iter = 1000000000;

	// results
	t->performance = -1.0;
	t->duration = -1.0;
	t->durations = NULL;
	t->predictions = NULL;

	return t;
}

/**
 * @brief Free the GenTask struct
 *
 * @details
 * Freeing the allocated memory of the GenTask means freeing _only_ the task 
 * itself. The datasets are not freed, as these are shared between all tasks.
 *
 * @param[in] 	t 	GenTask to be freed
 *
 */
void gensvm_free_task(struct GenTask *t)
{
	if (t->predictions != NULL)
		free(t->predictions);
	if (t->durations != NULL)
		free(t->durations);
	free(t);
	t = NULL;
}

/**
 * @brief Deepcopy a GenTask struct
 *
 * @details
 * Create a deep copy of a GenTask struct. Note that the datasets belonging to 
 * the tasks are not copied, only the pointers to the datasets.
 *
 * @param[in] 	t 	input GenTask struct to copy
 *
 * @return 		a deepcopy of the input GenTask
 */
struct GenTask *gensvm_copy_task(struct GenTask *t)
{
	struct GenTask *nt = gensvm_init_task();
	nt->weight_idx = t->weight_idx;
	nt->folds = t->folds;
	nt->ID = t->ID;
	nt->p = t->p;
	nt->kappa = t->kappa;
	nt->lambda = t->lambda;
	nt->epsilon = t->epsilon;
	nt->train_data = t->train_data;
	nt->test_data = t->test_data;
	nt->performance = t->performance;

	nt->kerneltype = t->kerneltype;
	nt->gamma = t->gamma;
	nt->coef = t->coef;
	nt->degree = t->degree;

	nt->max_iter = t->max_iter;

	return nt;
}

/**
 * @brief Copy parameters from GenTask to GenModel
 *
 * @details
 * A GenTask struct only contains the parameters of the GenModel to be estimated.
 * This function is used to copy these parameters.
 *
 * @param[in] 		task 	GenTask instance with parameters
 * @param[in,out] 	model 	GenModel to which the parameters are copied
 */
void gensvm_task_to_model(struct GenTask *task, struct GenModel *model)
{
	// copy basic model parameters
	model->weight_idx = task->weight_idx;
	model->epsilon = task->epsilon;
	model->p = task->p;
	model->kappa = task->kappa;
	model->lambda = task->lambda;

	// copy kernel parameters
	model->kerneltype = task->kerneltype;
	model->gamma = task->gamma;
	model->coef = task->coef;
	model->degree = task->degree;

	// copy other parameters
	model->max_iter = task->max_iter;
}
