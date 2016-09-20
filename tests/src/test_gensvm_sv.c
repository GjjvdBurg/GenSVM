/**
 * @file test_gensvm_sv.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Unit tests for gensvm_sv.c functions
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
