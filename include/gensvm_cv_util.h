/**
 * @file gensvm_cv_util.h
 * @author G.J.J. van den Burg
 * @date 2014-01-07
 * @brief Header file for gensvm_cv_util.c
 *
 * @details
 * Contains function declarations for functions needed for performing cross
 * validation on GenData structures.
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

#ifndef GENSVM_CV_UTIL_H
#define GENSVM_CV_UTIL_H

#include "gensvm_base.h"
#include "gensvm_rand.h"

void gensvm_make_cv_split(long N, long folds, long *cv_idx);
void gensvm_get_tt_split(struct GenData *full_data, struct GenData *train_data,
		struct GenData *test_data, long *cv_idx, long fold_idx);
void gensvm_get_tt_split_dense(struct GenData *full_data,
		struct GenData *train_data, struct GenData *test_data,
		long *cv_idx, long fold_idx);
void gensvm_get_tt_split_sparse(struct GenData *full_data,
		struct GenData *train_data, struct GenData *test_data,
		long *cv_idx, long fold_idx);

#endif
