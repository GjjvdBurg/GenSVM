/**
 * @file test_gensvm_base.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Unit tests for gensvm_base.c functions
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

	return NULL;
}

RUN_TESTS(all_tests);
