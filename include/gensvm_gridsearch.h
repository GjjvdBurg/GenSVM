/**
 * @file gensvm_gridsearch.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Header file for gensvm_gridsearch.c
 *
 * @details
 * The grid search for the optimal parameters is done through a queue.
 * This file contains struct definitions for this queue and a single
 * task in a queue, as well as a structure for the complete training
 * scheme. Function declarations are also included.
 *
 */

#ifndef GENSVM_GRIDSEARCH_H
#define GENSVM_GRIDSEARCH_H

// includes
#include "gensvm_cv_util.h"
#include "gensvm_init.h"
#include "gensvm_grid.h"
#include "gensvm_optimize.h"
#include "gensvm_pred.h"
#include "gensvm_queue.h"
#include "gensvm_timer.h"

// function declarations
void gensvm_fill_queue(struct GenGrid *grid, struct GenQueue *queue,
		struct GenData *train_data, struct GenData *test_data);
void consistency_repeats(struct GenQueue *q, long repeats, TrainType traintype);
void make_model_from_task(struct GenTask *task, struct GenModel *model);
void print_progress_string(struct GenTask *task, long N);
void start_training(struct GenQueue *q);
double gensvm_cross_validation(struct GenModel *model,
	       	struct GenData **train_folds, struct GenData **test_folds,
		int folds, long n_total);
#endif
