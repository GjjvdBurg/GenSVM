/**
 * @file gensvm_copy.c
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Function for copying a GenModel instance
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

#include "gensvm_copy.h"

/**
 * @brief Copy model parameters between two GenModel structs
 *
 * @details
 * The parameters copied are:
 *  - GenModel::weight_idx
 *  - GenModel::epsilon
 *  - GenModel::p
 *  - GenModel::kappa
 *  - GenModel::lambda
 *  - GenModel::kerneltype
 *  - GenModel::gamma
 *  - GenModel::coef
 *  - GenModel::degree
 *  - GenModel::max_iter
 *  - GenModel::seed
 *
 * @param[in] 		from 	GenModel to copy parameters from
 * @param[in,out] 	to 	GenModel to copy parameters to
 */
void gensvm_copy_model(struct GenModel *from, struct GenModel *to)
{
	to->weight_idx = from->weight_idx;
	to->epsilon = from->epsilon;
	to->p = from->p;
	to->kappa = from->kappa;
	to->lambda = from->lambda;

	to->kerneltype = from->kerneltype;
	to->gamma = from->gamma;
	to->coef = from->coef;
	to->degree = from->degree;

	to->max_iter = from->max_iter;
	to->seed = from->seed;
}
