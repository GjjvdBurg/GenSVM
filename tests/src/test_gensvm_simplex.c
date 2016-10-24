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

	mu_assert(fabs(matrix_get(model->U, 1, 0, 0) - -0.5) < 1e-14, 
			"U(0, 0) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 1, 1, 0) - 0.5) < 1e-14, 
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

	mu_assert(fabs(matrix_get(model->U, 3, 0, 0) - -0.5) < 1e-14,
				"U(0, 0) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 3, 1, 0) - 0.5) < 1e-14,
		"U(1, 0) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 3, 2, 0) - 0.0) < 1e-14,
		"U(2, 0) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 3, 3, 0) - 0.0) < 1e-14,
		"U(3, 0) incorrect.");

	mu_assert(fabs(matrix_get(model->U, 3, 0, 1) - -0.5/sqrt(3)) < 1e-14,
		  	"U(0, 1) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 3, 1, 1) - -0.5/sqrt(3)) < 1e-14,
		  	"U(1, 1) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 3, 2, 1) - 1.0/sqrt(3)) < 1e-14,
		  	"U(2, 1) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 3, 3, 1) - 0.0) < 1e-14,
		  	"U(3, 1) incorrect.");

	mu_assert(fabs(matrix_get(model->U, 3, 0, 2) - -1.0/sqrt(24)) < 1e-14,
		  	"U(0, 2) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 3, 1, 2) - -1.0/sqrt(24)) < 1e-14,
		  	"U(1, 2) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 3, 2, 2) - -1.0/sqrt(24)) < 1e-14,
		  	"U(2, 2) incorrect.");
	mu_assert(fabs(matrix_get(model->U, 3, 3, 2) - 3.0/sqrt(24)) < 1e-14,
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
	mu_assert(fabs(model->UU[0] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 0");
	mu_assert(fabs(model->UU[1] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 1");
	mu_assert(fabs(model->UU[2] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 2");
	mu_assert(fabs(model->UU[3] - -1.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 3");
	mu_assert(fabs(model->UU[4] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 4");
	mu_assert(fabs(model->UU[5] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 5");
	mu_assert(fabs(model->UU[6] - -0.5000000000000000) < 1e-14,
			"Incorrect value of model->UU at 6");
	mu_assert(fabs(model->UU[7] - -0.8660254037844388) < 1e-14,
			"Incorrect value of model->UU at 7");
	mu_assert(fabs(model->UU[8] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 8");
	mu_assert(fabs(model->UU[9] - -0.5000000000000000) < 1e-14,
			"Incorrect value of model->UU at 9");
	mu_assert(fabs(model->UU[10] - -0.2886751345948129) < 1e-14,
			"Incorrect value of model->UU at 10");
	mu_assert(fabs(model->UU[11] - -0.8164965809277261) < 1e-14,
			"Incorrect value of model->UU at 11");
	mu_assert(fabs(model->UU[12] - 1.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 12");
	mu_assert(fabs(model->UU[13] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 13");
	mu_assert(fabs(model->UU[14] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 14");
	mu_assert(fabs(model->UU[15] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 15");
	mu_assert(fabs(model->UU[16] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 16");
	mu_assert(fabs(model->UU[17] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 17");
	mu_assert(fabs(model->UU[18] - 0.5000000000000000) < 1e-14,
			"Incorrect value of model->UU at 18");
	mu_assert(fabs(model->UU[19] - -0.8660254037844388) < 1e-14,
			"Incorrect value of model->UU at 19");
	mu_assert(fabs(model->UU[20] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 20");
	mu_assert(fabs(model->UU[21] - 0.5000000000000000) < 1e-14,
			"Incorrect value of model->UU at 21");
	mu_assert(fabs(model->UU[22] - -0.2886751345948129) < 1e-14,
			"Incorrect value of model->UU at 22");
	mu_assert(fabs(model->UU[23] - -0.8164965809277261) < 1e-14,
			"Incorrect value of model->UU at 23");
	mu_assert(fabs(model->UU[24] - 0.5000000000000000) < 1e-14,
			"Incorrect value of model->UU at 24");
	mu_assert(fabs(model->UU[25] - 0.8660254037844388) < 1e-14,
			"Incorrect value of model->UU at 25");
	mu_assert(fabs(model->UU[26] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 26");
	mu_assert(fabs(model->UU[27] - -0.5000000000000000) < 1e-14,
			"Incorrect value of model->UU at 27");
	mu_assert(fabs(model->UU[28] - 0.8660254037844388) < 1e-14,
			"Incorrect value of model->UU at 28");
	mu_assert(fabs(model->UU[29] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 29");
	mu_assert(fabs(model->UU[30] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 30");
	mu_assert(fabs(model->UU[31] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 31");
	mu_assert(fabs(model->UU[32] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 32");
	mu_assert(fabs(model->UU[33] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 33");
	mu_assert(fabs(model->UU[34] - 0.5773502691896258) < 1e-14,
			"Incorrect value of model->UU at 34");
	mu_assert(fabs(model->UU[35] - -0.8164965809277261) < 1e-14,
			"Incorrect value of model->UU at 35");
	mu_assert(fabs(model->UU[36] - 0.5000000000000000) < 1e-14,
			"Incorrect value of model->UU at 36");
	mu_assert(fabs(model->UU[37] - 0.2886751345948129) < 1e-14,
			"Incorrect value of model->UU at 37");
	mu_assert(fabs(model->UU[38] - 0.8164965809277261) < 1e-14,
			"Incorrect value of model->UU at 38");
	mu_assert(fabs(model->UU[39] - -0.5000000000000000) < 1e-14,
			"Incorrect value of model->UU at 39");
	mu_assert(fabs(model->UU[40] - 0.2886751345948129) < 1e-14,
			"Incorrect value of model->UU at 40");
	mu_assert(fabs(model->UU[41] - 0.8164965809277261) < 1e-14,
			"Incorrect value of model->UU at 41");
	mu_assert(fabs(model->UU[42] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 42");
	mu_assert(fabs(model->UU[43] - -0.5773502691896258) < 1e-14,
			"Incorrect value of model->UU at 43");
	mu_assert(fabs(model->UU[44] - 0.8164965809277261) < 1e-14,
			"Incorrect value of model->UU at 44");
	mu_assert(fabs(model->UU[45] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 45");
	mu_assert(fabs(model->UU[46] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 46");
	mu_assert(fabs(model->UU[47] - 0.0000000000000000) < 1e-14,
			"Incorrect value of model->UU at 47");

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
