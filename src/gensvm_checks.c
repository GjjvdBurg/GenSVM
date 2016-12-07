/**
 * @file gensvm_checks.c
 * @author G.J.J. van den Burg
 * @date 2016-12-07
 * @brief Sanity checks used to ensure inputs are as expected
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

#include "gensvm_checks.h"

/**
 * @brief Check if the labels are contiguous on [1 .. K]
 *
 * @details
 * The GenSVM library currently requires that the labels that are supplied in 
 * a dataset are contigous on the interval [1 .. K] and have no gaps. This is 
 * required because the dimensionality of the problem is directly related to 
 * the maximum class label K. This function checks if the labels are indeed in 
 * the desired range.
 *
 * @param[in] 	data 	a GenData struct with the current data
 *
 * @return 		whether the labels are contiguous or not
 */
bool gensvm_check_outcome_contiguous(struct GenData *data)
{
	bool in_uniq, is_contiguous = true;
	long i, j, K = 1;
	long max_y = -1,
	     min_y = LONG_MAX;
	long *uniq_y = Calloc(long, K);
	uniq_y[0] = data->y[0];

	for (i=1; i<data->n; i++) {
		in_uniq = false;
		for (j=0; j<K; j++) {
			if (uniq_y[j] == data->y[i]) {
				in_uniq = true;
				break;
			}
		}

		if (!in_uniq) {
			uniq_y = Realloc(uniq_y, long, K+1);
			uniq_y[K++] = data->y[i];
		}

		max_y = maximum(max_y, data->y[i]);
		min_y = minimum(min_y, data->y[i]);
	}

	if (min_y < 1 || max_y > K) {
		is_contiguous = false;
	}

	free(uniq_y);

	return is_contiguous;
}
