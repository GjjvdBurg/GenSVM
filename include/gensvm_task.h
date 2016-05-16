/**
 * @file gensvm_task.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Struct for a single task in the queue
 *
 * @details
 * The grid search for the optimal parameters is done through a queue.
 * This file contains struct definitions for the tasks in the queue.
 * Initialization and free functions are also included.
 *
 */

#ifndef GENSVM_TASK_H
#define GENSVM_TASK_H

#include "gensvm_base.h"

/**
 * @brief A structure for a single task in the queue.
 *
 * @param folds 	number of folds in cross validation
 * @param ID 		numeric id of the task in the queue
 * @param weight_idx 	parameter for the GenModel
 * @param p 		parameter for the GenModel
 * @param kappa 	parameter for the GenModel
 * @param lambda 	parameter for the GenModel
 * @param epsilon 	parameter for the GenModel
 * @param kerneltype 	parameter for the GenModel
 * @param *kernelparam parameters for the GenModel
 * @param *train_data 	pointer to the training data
 * @param *test_data 	pointer to the test data (if any)
 * @param performance 	performance after cross validation
 */
struct GenTask {
	KernelType kerneltype;
	int weight_idx;
	long folds;
	long ID;
	double p;
	double kappa;
	double lambda;
	double epsilon;
	double *kernelparam;
	struct GenData *train_data;
	struct GenData *test_data;
	double performance;
};

struct GenTask *gensvm_init_task();
void gensvm_free_task(struct GenTask *task);

#endif
