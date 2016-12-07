/**
 * @file test_gensvm_checks.c
 * @author G.J.J. van den Burg
 * @date 2016-12-07
 * @brief Unit tests for gensvm_checks.c functions
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
#include "gensvm_checks.h"

char *test_check_outcome_contiguous_correct()
{
	struct GenData *data = gensvm_init_data();
	data->n = 10;
	data->y = Calloc(long, data->n);
	data->y[0] = 1;
	data->y[1] = 2;
	data->y[2] = 3;
	data->y[3] = 4;
	data->y[4] = 1;
	data->y[5] = 1;
	data->y[6] = 2;
	data->y[7] = 2;
	data->y[8] = 4;
	data->y[9] = 3;

	// start test code //
	mu_assert(gensvm_check_outcome_contiguous(data) == true, 
			"Incorrect check outcome for correct");
	// end test code //

	gensvm_free_data(data);

	return NULL;
}

char *test_check_outcome_contiguous_gap()
{
	struct GenData *data = gensvm_init_data();
	data->n = 10;
	data->y = Calloc(long, data->n);
	data->y[0] = 1;
	data->y[1] = 2;
	data->y[2] = 4;
	data->y[3] = 4;
	data->y[4] = 1;
	data->y[5] = 1;
	data->y[6] = 2;
	data->y[7] = 2;
	data->y[8] = 4;
	data->y[9] = 4;

	// start test code //
	mu_assert(gensvm_check_outcome_contiguous(data) == false, 
			"Incorrect check outcome for gap");
	// end test code //

	gensvm_free_data(data);

	return NULL;
}

char *test_check_outcome_contiguous_gaps()
{
	struct GenData *data = gensvm_init_data();
	data->n = 10;
	data->y = Calloc(long, data->n);
	data->y[0] = 1;
	data->y[1] = 6;
	data->y[2] = 4;
	data->y[3] = 4;
	data->y[4] = 1;
	data->y[5] = 1;
	data->y[6] = 6;
	data->y[7] = 6;
	data->y[8] = 4;
	data->y[9] = 4;

	// start test code //
	mu_assert(gensvm_check_outcome_contiguous(data) == false, 
			"Incorrect check outcome for gaps");
	// end test code //

	gensvm_free_data(data);

	return NULL;
}

char *test_check_outcome_contiguous_shift()
{
	struct GenData *data = gensvm_init_data();
	data->n = 10;
	data->y = Calloc(long, data->n);
	data->y[0] = 2;
	data->y[1] = 3;
	data->y[2] = 4;
	data->y[3] = 5;
	data->y[4] = 2;
	data->y[5] = 3;
	data->y[6] = 3;
	data->y[7] = 4;
	data->y[8] = 5;
	data->y[9] = 5;

	// start test code //
	mu_assert(gensvm_check_outcome_contiguous(data) == false, 
			"Incorrect check outcome for shift");
	// end test code //

	gensvm_free_data(data);

	return NULL;
}

char *all_tests()
{
	mu_suite_start();

	mu_run_test(test_check_outcome_contiguous_correct);
	mu_run_test(test_check_outcome_contiguous_gap);
	mu_run_test(test_check_outcome_contiguous_gaps);
	mu_run_test(test_check_outcome_contiguous_shift);

	return NULL;
}

RUN_TESTS(all_tests);
