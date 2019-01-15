/**
 * @file gensvm_optimize.h
 * @author G.J.J. van den Burg
 * @date 2013-08-01
 * @brief Header file for gensvm_optimize.c
 *
 * @details
 * Contains function declarations for functions used to train a single
 * GenModel.
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

#ifndef GENSVM_OPTIMIZE_H
#define GENSVM_OPTIMIZE_H

#include "gensvm_py_utils.h"
#include "gensvm_sv.h"
#include "gensvm_simplex.h"
#include "gensvm_timer.h"
#include "gensvm_update.h"
#include "gensvm_zv.h"

// function declarations
void gensvm_optimize(struct GenModel *model, struct GenData *data);
double gensvm_get_loss(struct GenModel *model, struct GenData *data, 
		struct GenWork *work);
void gensvm_calculate_errors(struct GenModel *model, struct GenData *data,
		double *ZV);
void gensvm_calculate_ZV_dense(struct GenModel *model, struct GenData *data,
		double *ZV);
void gensvm_calculate_ZV_sparse(struct GenModel *model, struct GenData *data,
		double *ZV);
void gensvm_calculate_huber(struct GenModel *model);
void gensvm_step_doubling(struct GenModel *model);

#endif
