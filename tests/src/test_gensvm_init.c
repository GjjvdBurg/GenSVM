/**
 * @file test_gensvm_init.c
 * @author G.J.J. van den Burg
 * @date 2016-08-01
 * @brief Unit tests for gensvm_init.c functions
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
#include "gensvm_init.h"

char *test_init_null_dense()
{
	int i;
	long n = 5,
	     m = 2,
	     K = 3;
	double value;
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	// start test code
	model->n = n;
	model->m = m;
	model->K = K;
	gensvm_allocate_model(model);

	data->n = n;
	data->m = m;
	data->K = K;
	data->RAW = Calloc(double, n*(m+1));
	for (i=0; i<n; i++) {
		matrix_set(data->RAW, data->n, m+1, i, 0, 1.0);
		matrix_set(data->RAW, data->n, m+1, i, 1, i);
		matrix_set(data->RAW, data->n, m+1, i, 2, -i);
	}
	data->Z = data->RAW;

	gensvm_init_V(NULL, model, data);

	// first row all ones
	value = matrix_get(model->V, m+1, K-1, 0, 0);
	mu_assert(value == 1.0, "Incorrect value at 0, 0");
	value = matrix_get(model->V, m+1, K-1, 0, 1);
	mu_assert(value == 1.0, "Incorrect value at 0, 1");

	// second row between -1 and 0.25
	value = matrix_get(model->V, m+1, K-1, 1, 0);
	mu_assert(value >= -1.0 && value <= 0.25, "Incorrect value at 1, 0");
	value = matrix_get(model->V, m+1, K-1, 1, 1);
	mu_assert(value >= -1.0 && value <= 0.25, "Incorrect value at 1, 1");

	// third row between -0.25 and 1
	value = matrix_get(model->V, m+1, K-1, 2, 0);
	mu_assert(value >= -0.25 && value <= 1.0, "Incorrect value at 2, 0");
	value = matrix_get(model->V, m+1, K-1, 2, 1);
	mu_assert(value >= -0.25 && value <= 1.0, "Incorrect value at 2, 1");

	// end test code

	gensvm_free_model(model);
	gensvm_free_data(data);

	return NULL;
}

char *test_init_null_sparse()
{
	int i;
	long n = 5,
	     m = 2,
	     K = 3;
	double value;
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	// start test code
	model->n = n;
	model->m = m;
	model->K = K;
	gensvm_allocate_model(model);

	data->n = n;
	data->m = m;
	data->K = K;
	data->RAW = Calloc(double, n*(m+1));
	for (i=0; i<n; i++) {
		matrix_set(data->RAW, n, m+1, i, 0, 1.0);
		matrix_set(data->RAW, n, m+1, i, 1, i);
		matrix_set(data->RAW, n, m+1, i, 2, -i);
	}

	data->Z = data->RAW;
	data->spZ = gensvm_dense_to_sparse(data->Z, n, m+1);
	free(data->RAW);
	data->RAW = NULL;
	data->Z = NULL;

	gensvm_init_V(NULL, model, data);

	// first row all ones
	value = matrix_get(model->V, m+1, K-1, 0, 0);
	mu_assert(value == 1.0, "Incorrect value at 0, 0");
	value = matrix_get(model->V, m+1, K-1, 0, 1);
	mu_assert(value == 1.0, "Incorrect value at 0, 1");

	// second row between -1 and 0.25
	value = matrix_get(model->V, m+1, K-1, 1, 0);
	mu_assert(value >= -1.0 && value <= 0.25, "Incorrect value at 1, 0");
	value = matrix_get(model->V, m+1, K-1, 1, 1);
	mu_assert(value >= -1.0 && value <= 0.25, "Incorrect value at 1, 1");

	// third row between -0.25 and 1
	value = matrix_get(model->V, m+1, K-1, 2, 0);
	mu_assert(value >= -0.25 && value <= 1.0, "Incorrect value at 2, 0");
	value = matrix_get(model->V, m+1, K-1, 2, 1);
	mu_assert(value >= -0.25 && value <= 1.0, "Incorrect value at 2, 1");

	// end test code

	gensvm_free_model(model);
	gensvm_free_data(data);

	return NULL;
}

char *test_init_seed()
{
	long n = 7,
	     m = 5,
	     K = 3;
	struct GenModel *model = gensvm_init_model();
	struct GenModel *seed = gensvm_init_model();
	struct GenData *data = gensvm_init_data();

	model->n = n;
	model->m = m;
	model->K = K;
	seed->n = n;
	seed->m = m;
	seed->K = K;
	data->n = n;
	data->m = m;
	data->K = K;
	gensvm_allocate_model(model);
	gensvm_allocate_model(seed);

	// start test code
	matrix_set(seed->V, seed->m+1, seed->K-1, 0, 0, 123.0);
	matrix_set(seed->V, seed->m+1, seed->K-1, 1, 1, 321.0);
	matrix_set(seed->V, seed->m+1, seed->K-1, 2, 0, 111.0);
	matrix_set(seed->V, seed->m+1, seed->K-1, 5, 0, 222.0);
	matrix_set(seed->V, seed->m+1, seed->K-1, 3, 1, 333.0);

	gensvm_init_V(seed, model, data);

	mu_assert(matrix_get(model->V, model->m+1, model->K-1, 0, 0) == 123.0,
			"Incorrect V value at 0, 0");
	mu_assert(matrix_get(model->V, model->m+1, model->K-1, 0, 1) == 0.0,
			"Incorrect V value at 0, 1");
	mu_assert(matrix_get(model->V, model->m+1, model->K-1, 1, 0) == 0.0,
			"Incorrect V value at 1, 0");
	mu_assert(matrix_get(model->V, model->m+1, model->K-1, 1, 1) == 321.0,
			"Incorrect V value at 1, 1");
	mu_assert(matrix_get(model->V, model->m+1, model->K-1, 2, 0) == 111.0,
			"Incorrect V value at 2, 0");
	mu_assert(matrix_get(model->V, model->m+1, model->K-1, 2, 1) == 0.0,
			"Incorrect V value at 2, 1");
	mu_assert(matrix_get(model->V, model->m+1, model->K-1, 3, 0) == 0.0,
			"Incorrect V value at 3, 0");
	mu_assert(matrix_get(model->V, model->m+1, model->K-1, 3, 1) == 333.0,
			"Incorrect V value at 3, 1");
	mu_assert(matrix_get(model->V, model->m+1, model->K-1, 4, 0) == 0.0,
			"Incorrect V value at 4, 0");
	mu_assert(matrix_get(model->V, model->m+1, model->K-1, 4, 1) == 0.0,
			"Incorrect V value at 4, 1");
	mu_assert(matrix_get(model->V, model->m+1, model->K-1, 5, 0) == 222.0,
			"Incorrect V value at 5, 0");
	mu_assert(matrix_get(model->V, model->m+1, model->K-1, 5, 1) == 0.0,
			"Incorrect V value at 5, 1");
	// end test code

	gensvm_free_model(model);
	gensvm_free_model(seed);
	gensvm_free_data(data);

	return NULL;
}

char *test_init_weights_1()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = NULL;
	model->n = 7;
	model->m = 5;
	model->K = 3;
	model->weight_idx = 1;
	gensvm_allocate_model(model);

	// start test code
	int i;
	gensvm_initialize_weights(data, model);
	for (i=0; i<model->n; i++) {
		mu_assert(model->rho[i] == 1.0, "incorrect weight in rho");
	}
	// end test code

	gensvm_free_model(model);
	gensvm_free_data(data);

	return NULL;
}

char *test_init_weights_2()
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();
	model->n = 8;
	model->m = 5;
	model->K = 3;
	model->weight_idx = 2;
	gensvm_allocate_model(model);

	data->y = Calloc(long, model->n);

	data->y[0] = 1;
	data->y[1] = 1;
	data->y[2] = 1;
	data->y[3] = 1;
	data->y[4] = 2;
	data->y[5] = 2;
	data->y[6] = 2;
	data->y[7] = 3;

	// start test code
	gensvm_initialize_weights(data, model);
	int i;
	for (i=0; i<4; i++) {
		mu_assert(model->rho[i] == 8.0/(4.0 * 3.0),
				"Incorrect weight for class 1");
	}
	for (i=0; i<3; i++) {
		mu_assert(model->rho[4+i] == 8.0/(3.0 * 3.0),
				"Incorrect weight for class 2");
	}
	mu_assert(model->rho[7] == 8.0/(1.0 * 3.0),
			"Incorrect weight for class 3");

	// end test code

	gensvm_free_model(model);
	gensvm_free_data(data);

	return NULL;
}

char *test_init_weights_wrong()
{
	mu_test_missing();
	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_init_null_dense);
	mu_run_test(test_init_null_sparse);
	mu_run_test(test_init_seed);
	mu_run_test(test_init_weights_1);
	mu_run_test(test_init_weights_2);
	mu_run_test(test_init_weights_wrong);

	return NULL;
}

RUN_TESTS(all_tests);
