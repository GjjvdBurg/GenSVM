/**
 * @file test_gensvm_sv.c
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Unit tests for gensvm_sv.c functions
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

#include "minunit.h"
#include "gensvm_sv.h"

char *test_sv()
{
	struct GenModel *model = gensvm_init_model();

	model->n = 5;
	model->m = 3;
	model->K = 3;
	gensvm_allocate_model(model);

	// for a support vector we need less than 2 elements per row larger 
	// than 1

	// this is an sv
	matrix_set(model->Q, model->K, 0, 0, 1.1);
	matrix_set(model->Q, model->K, 0, 1, 0.0);
	matrix_set(model->Q, model->K, 0, 2, 1.0);

	// this is an sv
	matrix_set(model->Q, model->K, 1, 0, 0.5);
	matrix_set(model->Q, model->K, 1, 1, 0.5);
	matrix_set(model->Q, model->K, 1, 2, 0.5);

	// this is an sv
	matrix_set(model->Q, model->K, 2, 0, -0.5);
	matrix_set(model->Q, model->K, 2, 1, 0.5);
	matrix_set(model->Q, model->K, 2, 2, -0.5);

	// this is not an sv
	matrix_set(model->Q, model->K, 3, 0, 1.5);
	matrix_set(model->Q, model->K, 3, 1, 1.5);
	matrix_set(model->Q, model->K, 3, 2, 0.5);

	// this is not an sv
	matrix_set(model->Q, model->K, 4, 0, 2.0);
	matrix_set(model->Q, model->K, 4, 1, 2.0);
	matrix_set(model->Q, model->K, 4, 2, 2.0);

	mu_assert(gensvm_num_sv(model) == 3, "number of svs incorrect");

	gensvm_free_model(model);
	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_sv);

	return NULL;
}

RUN_TESTS(all_tests);
