/**
 * @file test_gensvm_base.c
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Unit tests for gensvm_base.c functions
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
#include "gensvm_base.h"

char *test_init_free_model_1()
{
	struct GenModel *model = gensvm_init_model();
	gensvm_free_model(model);
	return NULL;
}

char *test_init_free_model_2()
{
	struct GenModel *model = gensvm_init_model();
	gensvm_free_model(model);
	gensvm_free_model(NULL);
	return NULL;
}

char *test_allocate_free_model()
{
	struct GenModel *model = gensvm_init_model();
	model->n = 3;
	model->m = 4;
	model->K = 5;
	gensvm_allocate_model(model);
	gensvm_free_model(model);
	return NULL;
}

char *test_reallocate_free_model()
{
	struct GenModel *model = gensvm_init_model();
	model->n = 3;
	model->m = 4;
	model->K = 5;
	gensvm_allocate_model(model);
	gensvm_reallocate_model(model, 3, 4);
	model->n = 4;
	gensvm_reallocate_model(model, 4, 4);
	gensvm_reallocate_model(model, 4, 5);
	gensvm_reallocate_model(model, 3, 4);

	gensvm_free_model(model);
	return NULL;
}

char *test_init_free_data_1()
{
	struct GenData *data = gensvm_init_data();
	gensvm_free_data(data);
	return NULL;
}

char *test_init_free_data_2()
{
	struct GenData *data = NULL;
	gensvm_free_data(data);
	return NULL;
}

char *test_init_free_data_3()
{
	struct GenData *data = gensvm_init_data();
	data->Z = Calloc(double, 3);
	gensvm_free_data(data);
	return NULL;
}

char *test_init_free_work()
{
	struct GenModel *model = gensvm_init_model();
	model->n = 10;
	model->m = 4;
	model->K = 3;

	struct GenWork *work = gensvm_init_work(model);

	mu_assert(model->n == work->n, "n variable copied incorrectly");
	mu_assert(model->m == work->m, "m variable copied incorrectly");
	mu_assert(model->K == work->K, "K variable copied incorrectly");

	mu_assert(work->LZ != NULL, "LZ variable is NULL");
	mu_assert(work->ZB != NULL, "ZB variable is NULL");
	mu_assert(work->ZBc != NULL, "ZBc variable is NULL");
	mu_assert(work->ZAZ != NULL, "ZAZ variable is NULL");
	mu_assert(work->ZV != NULL, "ZV variable is NULL");
	mu_assert(work->beta != NULL, "beta variable is NULL");

	gensvm_free_model(model);
	gensvm_free_work(work);

	return NULL;
}

void fill_with_noise(double *ptr, int elem)
{
	int i;
	for (i=0; i<elem; i++) {
		ptr[i] = ((double) rand())/((double) RAND_MAX);
	}
}

bool all_elements_zero(double *ptr, int elem)
{
	int i;
	for (i=0; i<elem; i++) {
		if (ptr[i] != 0)
			return false;
	}
	return true;
}

char *test_reset_work()
{
	struct GenModel *model = gensvm_init_model();
	int n = 10,
	    m = 4,
	    K = 3;

	model->n = n;
	model->m = m;
	model->K = K;

	struct GenWork *work = gensvm_init_work(model);

	// start test code //
	fill_with_noise(work->LZ, n*(m+1));
	fill_with_noise(work->ZB, (m+1)*(K-1));
	fill_with_noise(work->ZBc, (m+1)*(K-1));
	fill_with_noise(work->ZAZ, (m+1)*(m+1));
	fill_with_noise(work->ZV, n*(K-1));
	fill_with_noise(work->beta, K-1);

	gensvm_reset_work(work);

	mu_assert(all_elements_zero(work->LZ, n*(m+1)), 
			"Not all elements of LZ are zero");
	mu_assert(all_elements_zero(work->ZB, (m+1)*(K-1)),
			"Not all elements of ZB are zero");
	mu_assert(all_elements_zero(work->ZBc, (m+1)*(K-1)),
				"Not all elements of ZBc are zero");
	mu_assert(all_elements_zero(work->ZAZ, (m+1)*(m+1)),
				"Not all elements of ZAZ are zero");
	mu_assert(all_elements_zero(work->ZV, n*(K-1)),
			"Not all elements of ZV are zero");
	mu_assert(all_elements_zero(work->beta, K-1),
			"Not all elements of beta are zero");

	// end test code //

	gensvm_free_model(model);
	gensvm_free_work(work);

	return NULL;
}


char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_init_free_model_1);
	mu_run_test(test_init_free_model_2);
	mu_run_test(test_allocate_free_model);
	mu_run_test(test_reallocate_free_model);

	mu_run_test(test_init_free_data_1);
	mu_run_test(test_init_free_data_2);
	mu_run_test(test_init_free_data_3);

	mu_run_test(test_init_free_work);
	mu_run_test(test_reset_work);

	return NULL;
}

RUN_TESTS(all_tests);
