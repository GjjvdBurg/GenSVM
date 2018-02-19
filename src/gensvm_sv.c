/**
 * @file gensvm_sv.c
 * @author G.J.J. van den Burg
 * @date 2014-05-01
 * @brief Calculate the number of support vectors
 *
 * @details
 * The function in this file can be used to calculate the number of support
 * vectors are left in a model.
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

#include "gensvm_sv.h"

/**
 * @brief Calculate the number of support vectors in a model
 *
 * @details
 * If an object is correctly classified, the number of classes for which the
 * error q is larger than 1, is K-1 (i.e., there is no error w.r.t. any of the
 * other classes).  All objects for which this is not the case are thus support
 * vectors.
 *
 * @param[in] 	model 	GenModel with solution and up-to-date Q matrix
 * @return 		number of support vectors with this solution
 *
 */
long gensvm_num_sv(struct GenModel *model)
{
	long i, j, num_correct, num_sv = 0;
	double value;

	for (i=0; i<model->n; i++) {
		num_correct = 0;
		for (j=0; j<model->K; j++) {
			value = matrix_get(model->Q, model->n, model->K, i, j);
			num_correct += (value > 1);
		}
		num_sv += (num_correct < model->K - 1);
	}

	return num_sv;
}
