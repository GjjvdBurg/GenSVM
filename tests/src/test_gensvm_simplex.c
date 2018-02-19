/**
 * @file test_gensvm_simplex.c
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Unit tests for gensvm_simplex.c functions
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
#include "gensvm_simplex.h"

char *test_simplex_1()
{
	struct GenModel *model = gensvm_init_model();
	model->U = Calloc(double, 2*1);
	model->K = 2;

	gensvm_simplex(model);

	mu_assert(fabs(matrix_get(model->U, 2, 1, 0, 0) - -0.5) < 1e-14,
			"U(0, 0) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 2, 1, 1, 0) - 0.5) < 1e-14,
			"U(1, 0) incorrect.");

	gensvm_free_model(model);

	return NULL;
}

char *test_simplex_2()
{
	struct GenModel *model = gensvm_init_model();
	model->U = Calloc(double, 4*3);
	model->K = 4;

	gensvm_simplex(model);

	mu_assert(fabs(matrix_get(model->U, 4, 3, 0, 0) - -0.5) < 1e-14,
				"U(0, 0) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 4, 3, 1, 0) - 0.5) < 1e-14,
		"U(1, 0) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 4, 3, 2, 0) - 0.0) < 1e-14,
		"U(2, 0) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 4, 3, 3, 0) - 0.0) < 1e-14,
		"U(3, 0) incorrect.");

	mu_assert(fabs(matrix_get(model->U, 4, 3, 0, 1) - -0.5/sqrt(3)) < 1e-14,
			"U(0, 1) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 4, 3, 1, 1) - -0.5/sqrt(3)) < 1e-14,
			"U(1, 1) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 4, 3, 2, 1) - 1.0/sqrt(3)) < 1e-14,
			"U(2, 1) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 4, 3, 3, 1) - 0.0) < 1e-14,
			"U(3, 1) incorrect.");

	mu_assert(fabs(matrix_get(model->U, 4, 3, 0, 2) - -1.0/sqrt(24)) < 1e-14,
			"U(0, 2) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 4, 3, 1, 2) - -1.0/sqrt(24)) < 1e-14,
			"U(1, 2) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 4, 3, 2, 2) - -1.0/sqrt(24)) < 1e-14,
			"U(2, 2) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 4, 3, 3, 2) - 3.0/sqrt(24)) < 1e-14,
			"U(3, 2) incorrect.");

	gensvm_free_model(model);

	return NULL;
}

char *test_gensvm_simplex_diff()
{
	struct GenData *data = gensvm_init_data();
	struct GenModel *model = gensvm_init_model();

	int n = 8,
	    m = 3,
	    K = 4;
	model->n = n;
	model->m = m;
	model->K = K;
	data->n = n;
	data->m = m;
	data->K = K;

	gensvm_allocate_model(model);
	gensvm_simplex(model);

	// start test code //
	gensvm_simplex_diff(model);

	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 0, 0) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 0,0");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 0, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 0,1");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 0, 2) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 0,2");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 1, 0) -
				-1.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 1,0");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 1, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 1,1");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 1, 2) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 1,2");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 2, 0) -
				-0.5000000000000000) < 1e-14,
			"Incorrect value of model->UU at 2,0");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 2, 1) -
				-0.8660254037844388) < 1e-14,
			"Incorrect value of model->UU at 2,1");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 2, 2) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 2,2");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 3, 0) -
				-0.5000000000000000) < 1e-14,
			"Incorrect value of model->UU at 3,0");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 3, 1) -
				-0.2886751345948129) < 1e-14,
			"Incorrect value of model->UU at 3,1");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 3, 2) -
				-0.8164965809277261) < 1e-14,
			"Incorrect value of model->UU at 3,2");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 4, 0) -
				1.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 4,0");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 4, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 4,1");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 4, 2) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 4,2");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 5, 0) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 5,0");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 5, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 5,1");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 5, 2) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 5,2");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 6, 0) -
				0.5000000000000000) < 1e-14,
			"Incorrect value of model->UU at 6,0");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 6, 1) -
				-0.8660254037844388) < 1e-14,
			"Incorrect value of model->UU at 6,1");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 6, 2) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 6,2");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 7, 0) -
				0.5000000000000000) < 1e-14,
			"Incorrect value of model->UU at 7,0");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 7, 1) -
				-0.2886751345948129) < 1e-14,
			"Incorrect value of model->UU at 7,1");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 7, 2) -
				-0.8164965809277261) < 1e-14,
			"Incorrect value of model->UU at 7,2");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 8, 0) -
				0.5000000000000000) < 1e-14,
			"Incorrect value of model->UU at 8,0");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 8, 1) -
				0.8660254037844388) < 1e-14,
			"Incorrect value of model->UU at 8,1");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 8, 2) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 8,2");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 9, 0) -
				-0.5000000000000000) < 1e-14,
			"Incorrect value of model->UU at 9,0");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 9, 1) -
				0.8660254037844388) < 1e-14,
			"Incorrect value of model->UU at 9,1");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 9, 2) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 9,2");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 10, 0) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 10,0");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 10, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 10,1");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 10, 2) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 10,2");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 11, 0) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 11,0");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 11, 1) -
				0.5773502691896258) < 1e-14,
			"Incorrect value of model->UU at 11,1");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 11, 2) -
				-0.8164965809277261) < 1e-14,
			"Incorrect value of model->UU at 11,2");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 12, 0) -
				0.5000000000000000) < 1e-14,
			"Incorrect value of model->UU at 12,0");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 12, 1) -
				0.2886751345948129) < 1e-14,
			"Incorrect value of model->UU at 12,1");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 12, 2) -
				0.8164965809277261) < 1e-14,
			"Incorrect value of model->UU at 12,2");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 13, 0) -
				-0.5000000000000000) < 1e-14,
			"Incorrect value of model->UU at 13,0");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 13, 1) -
				0.2886751345948129) < 1e-14,
			"Incorrect value of model->UU at 13,1");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 13, 2) -
				0.8164965809277261) < 1e-14,
			"Incorrect value of model->UU at 13,2");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 14, 0) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 14,0");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 14, 1) -
				-0.5773502691896258) < 1e-14,
			"Incorrect value of model->UU at 14,1");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 14, 2) -
				0.8164965809277261) < 1e-14,
			"Incorrect value of model->UU at 14,2");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 15, 0) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 15,0");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 15, 1) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 15,1");
	mu_assert(fabs(matrix_get(model->UU, K*K, K-1, 15, 2) -
				0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 15,2");
	// end test code //

	gensvm_free_model(model);
	gensvm_free_data(data);

	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_simplex_1);
	mu_run_test(test_simplex_2);
	mu_run_test(test_gensvm_simplex_diff);

	return NULL;
}

RUN_TESTS(all_tests);
