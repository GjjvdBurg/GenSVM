/**
 * @file gensvm_consistency.h
 * @author G.J.J. van den Burg
 * @date 2016-10-24
 * @brief Header file for gensvm_consistency.c
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

#ifndef GENSVM_CONSISTENCY_H
#define GENSVM_CONSISTENCY_H

// includes
#include "gensvm_queue.h"
#include "gensvm_print.h"
#include "gensvm_cv_util.h"
#include "gensvm_cross_validation.h"
#include "gensvm_timer.h"

// function declarations
struct GenQueue *gensvm_top_queue(struct GenQueue *q, double percentile);
int gensvm_dsort(const void *elem1, const void *elem2);
int gensvm_consistency_repeats(struct GenQueue *q, long repeats,
		double percentile);
double gensvm_percentile(double *values, long N, double p);

#endif
