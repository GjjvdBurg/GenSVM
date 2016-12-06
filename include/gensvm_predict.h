/**
 * @file gensvm_predict.h
 * @author G.J.J. van den Burg
 * @date 2013-08-01
 * @brief Header file for gensvm_predict.c
 *
 * @details
 * Contains function declarations for prediction functions.
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

#ifndef GENSVM_PREDICT_H
#define GENSVM_PREDICT_H

// includes
#include "gensvm_kernel.h"
#include "gensvm_simplex.h"
#include "gensvm_zv.h"

// function declarations
void gensvm_predict_labels(struct GenData *testdata,
	       	struct GenModel *model, long *predy);
double gensvm_prediction_perf(struct GenData *data, long *perdy);

#endif
