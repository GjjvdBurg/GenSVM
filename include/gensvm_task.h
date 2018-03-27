/**
 * @file gensvm_task.h
 * @author G.J.J. van den Burg
 * @date 2013-08-01
 * @brief Header file for gensvm_task.c
 *
 * @details
 * The grid search for the optimal parameters is done through a queue.
 * This file contains struct definitions for the tasks in the queue.
 * Initialization and free functions are also included.
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

#ifndef GENSVM_TASK_H
#define GENSVM_TASK_H

#include "gensvm_base.h"

/**
 * @brief A structure for a single task in the queue.
 *
 * @param ID 		numeric id of the task in the queue
 * @param folds 	number of folds in cross validation
 * @param train_data 	pointer to the training data
 * @param test_data 	pointer to the test data (if any)
 *
 * @param kerneltype 	parameter for the GenModel
 * @param weight_idx 	parameter for the GenModel
 * @param p 		parameter for the GenModel
 * @param kappa 	parameter for the GenModel
 * @param lambda 	parameter for the GenModel
 * @param epsilon 	parameter for the GenModel
 * @param gamma 	parameter for the GenModel
 * @param coef 		parameter for the GenModel
 * @param degree 	parameter for the GenModel
 *
 * @param performance 	accuracy on training set during cross validation
 * @param duration 	training time in seconds
 * @param durations 	array of training time per fold in seconds
 * @param predictions 	array of predictions on the training set
 *
 */
struct GenTask {
	long ID;
	///< numeric id of the task in the queue
	long folds;
	///< number of folds in cross validation
	struct GenData *train_data;
	///< pointer to the training data
	struct GenData *test_data;
	///< pointer to the test data (if any)
	KernelType kerneltype;
	///< kerneltype parameter for the GenModel
	int weight_idx;
	///< weight_idx parameter for the GenModel
	double p;
	///< p parameter for the GenModel
	double kappa;
	///< kappa parameter for the GenModel
	double lambda;
	///< lambda parameter for the GenModel
	double epsilon;
	///< epsilon parameter for the GenModel
	double gamma;
	///< gamma parameter for the GenModel
	double coef;
	///< coef parameter for the GenModel
	double degree;
	///< degree parameter for the GenModel
	long max_iter;
	///< maximum number of iterations of the algorithm
	double performance;
	///< performance after cross validation
	double duration;
	///< training time in seconds
	double *durations;
	///< array of training time per fold in seconds
	long *predictions;
	///< array of CV predictions on the training data
};

struct GenTask *gensvm_init_task(void);
struct GenTask *gensvm_copy_task(struct GenTask *t);
void gensvm_free_task(struct GenTask *t);
void gensvm_task_to_model(struct GenTask *task, struct GenModel *model);

#endif
