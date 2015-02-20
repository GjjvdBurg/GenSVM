/**
 * @file gensvm_train_dataset.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Structs and functions necessary for the grid search
 *
 * @details
 * The grid search for the optimal parameters is done through a queue.
 * This file contains struct definitions for this queue and a single
 * task in a queue, as well as a structure for the complete training
 * scheme. Function declarations are also included.
 *
 */

#ifndef GENSVM_TRAIN_DATASET_H
#define GENSVM_TRAIN_DATASET_H

#include "globals.h"
#include "types.h"

// forward declarations
struct GenData;
struct GenModel;

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
struct Task {
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

/**
 * @brief Simple task queue.
 *
 * This struct is basically just an array of pointers to Task instances,
 * with a length and an index of the current task.
 *
 * @param **tasks 	array of pointers to Task structs
 * @param N 		size of task array
 * @param i 		index used for keeping track of the queue
 */
struct Queue {
	struct Task **tasks;
	long N;
	long i;
};

/**
 * @brief Structure for describing the entire grid search
 *
 * @param traintype 		type of training to use
 * @param kerneltype 		type of kernel to use throughout training
 * @param repeats 		number of repeats to be done after the grid 
 * 				search to find the parameter set with the 
 * 				most consistent high performance
 * @param folds 		number of folds in cross validation
 * @param Np 			size of the array of p values
 * @param Nl 			size of the array of lambda values
 * @param Nk 			size of the array of kappa values
 * @param Ne 			size of the array of epsilon values
 * @param Nw 			size of the array of weight_idx values
 * @param Ng 			size of the array of gamma values
 * @param Nc 			size of the array of coef values
 * @param Nd 			size of the array of degree values
 * @param *weight_idxs 		array of weight_idxs
 * @param *ps 			array of p values 
 * @param *lambdas 		array of lambda values
 * @param *kappas 		array of kappa values
 * @param *epsilons 		array of epsilon values
 * @param *gammas 		array of gamma values
 * @param *coefs 		array of coef values
 * @param *degrees 		array of degree values
 * @param *train_data_file 	filename of train data file
 * @param *test_data_file 	filename of test data file
 *
 */
struct Training {
	TrainType traintype;
	KernelType kerneltype;
	long repeats;
	long folds;
	long Np;
	long Nl;
	long Nk;
	long Ne;
	long Nw;
	long Ng;
	long Nc;
	long Nd;
	int *weight_idxs;
	double *ps;
	double *lambdas;
	double *kappas;
	double *epsilons;
	double *gammas;
	double *coefs;
	double *degrees;
	char *train_data_file;
	char *test_data_file;
};

void make_queue(struct Training *training, struct Queue *queue,
		struct GenData *train_data, struct GenData *test_data);

struct Task *get_next_task(struct Queue *q);
void start_training_tt(struct Queue *q);
void start_training_cv(struct Queue *q);
void free_queue(struct Queue *q);

void consistency_repeats(struct Queue *q, long repeats, TrainType traintype);

double cross_validation(struct GenModel *model, struct GenData *data,
	       	long folds);

void make_model_from_task(struct Task *task, struct GenModel *model);
void copy_model(struct GenModel *from, struct GenModel *to);

void print_progress_string(struct Task *task, long N);

// new
void start_training(struct Queue *q);
double gensvm_cross_validation(struct GenModel *model,
	       	struct GenData **train_folds, struct GenData **test_folds,
		int folds, long n_total);
#endif
