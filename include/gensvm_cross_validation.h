/**
 * @file gensvm_cross_validation.h
 * @author G.J.J. van den Burg
 * @date 2016-10-24
 * @brief Header file for gensvm_cross_validation.c
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

#ifndef GENSVM_CROSS_VALIDATION_H
#define GENSVM_CROSS_VALIDATION_H

// includes
#include "gensvm_base.h"
#include "gensvm_init.h"
#include "gensvm_optimize.h"
#include "gensvm_predict.h"
#include "gensvm_timer.h"

// function declarations
double gensvm_cross_validation(struct GenModel *model,
		struct GenData **train_folds, struct GenData **test_folds,
		long folds, long n_total);
void gensvm_cross_validation_store(struct GenModel *model, 
		struct GenData **train_folds, struct GenData **test_folds, 
		long folds, long n_total, long *cv_idx, long *predictions, 
		double *durations, int verbosity);

#endif
