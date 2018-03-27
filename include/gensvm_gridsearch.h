/**
 * @file gensvm_gridsearch.h
 * @author G.J.J. van den Burg
 * @date 2013-08-01
 * @brief Header file for gensvm_gridsearch.c
 *
 * @details
 * The grid search for the optimal parameters is done through a queue.
 * This file contains struct definitions for this queue and a single
 * task in a queue, as well as a structure for the complete training
 * scheme. Function declarations are also included.
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

#ifndef GENSVM_GRIDSEARCH_H
#define GENSVM_GRIDSEARCH_H

// includes
#include "gensvm_cross_validation.h"
#include "gensvm_cv_util.h"
#include "gensvm_grid.h"
#include "gensvm_queue.h"
#include "gensvm_timer.h"

// function declarations
void gensvm_fill_queue(struct GenGrid *grid, struct GenQueue *queue,
		struct GenData *train_data, struct GenData *test_data);
bool gensvm_kernel_changed(struct GenTask *newtask, struct GenTask *oldtask);
void gensvm_kernel_folds(long folds, struct GenModel *model,
		struct GenData **train_folds, struct GenData **test_folds);
void gensvm_gridsearch_progress(struct GenTask *task, long N, double perf,
		double duration, double current_max, bool show_perf);
double gensvm_train_queue(struct GenQueue *q, long *cv_idx,
		bool store_predictions, int verbosity);

#endif
