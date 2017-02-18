/**
 * @file test_gensvm_task.c
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Unit tests for gensvm_task.c functions
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
#include "gensvm_task.h"

char *test_init_free_task()
{
	struct GenTask *task = gensvm_init_task();
	gensvm_free_task(task);
	return NULL;
}

char *test_task_to_model_linear()
{
	struct GenTask *task = gensvm_init_task();
	struct GenModel *model = gensvm_init_model();

	// start test code //
	task->weight_idx = 2;
	task->p = 1.3;
	task->kappa = 0.1;
	task->lambda = 1.4;
	task->epsilon = 5e-3;
	task->kerneltype = K_LINEAR;
	task->max_iter = 100;

	gensvm_task_to_model(task, model);

	mu_assert(model->weight_idx = 2, "Incorrect model weight_idx");
	mu_assert(model->p == 1.3, "Incorrect model p");
	mu_assert(model->kappa == 0.1, "Incorrect model kappa");
	mu_assert(model->lambda == 1.4, "Incorrect model lambda");
	mu_assert(model->epsilon == 5e-3, "Incorrect model epsilon");
	mu_assert(model->kerneltype == K_LINEAR, "Incorrect model kerneltype");
	mu_assert(model->max_iter == 100, "Incorrect model max_iter");
	// end test code //

	gensvm_free_model(model);
	gensvm_free_task(task);
	return NULL;
}

char *test_task_to_model_rbf()
{
	struct GenTask *task = gensvm_init_task();
	struct GenModel *model = gensvm_init_model();

	// start test code //
	task->weight_idx = 2;
	task->p = 1.3;
	task->kappa = 0.1;
	task->lambda = 1.4;
	task->epsilon = 5e-3;
	task->kerneltype = K_RBF;
	task->gamma = 3.1;

	gensvm_task_to_model(task, model);

	mu_assert(model->weight_idx = 2, "Incorrect model weight_idx");
	mu_assert(model->p == 1.3, "Incorrect model p");
	mu_assert(model->kappa == 0.1, "Incorrect model kappa");
	mu_assert(model->lambda == 1.4, "Incorrect model lambda");
	mu_assert(model->epsilon == 5e-3, "Incorrect model epsilon");
	mu_assert(model->kerneltype == K_RBF, "Incorrect model kerneltype");
	mu_assert(model->gamma == 3.1, "Incorrect model gamma");
	// end test code //

	gensvm_free_model(model);
	gensvm_free_task(task);
	return NULL;
}

char *test_task_to_model_poly()
{
	struct GenTask *task = gensvm_init_task();
	struct GenModel *model = gensvm_init_model();

	// start test code //
	task->weight_idx = 2;
	task->p = 1.3;
	task->kappa = 0.1;
	task->lambda = 1.4;
	task->epsilon = 5e-3;
	task->kerneltype = K_POLY;
	task->gamma = 3.1;
	task->coef = 2.1;
	task->degree = 1.1;

	gensvm_task_to_model(task, model);

	mu_assert(model->weight_idx = 2, "Incorrect model weight_idx");
	mu_assert(model->p == 1.3, "Incorrect model p");
	mu_assert(model->kappa == 0.1, "Incorrect model kappa");
	mu_assert(model->lambda == 1.4, "Incorrect model lambda");
	mu_assert(model->epsilon == 5e-3, "Incorrect model epsilon");
	mu_assert(model->kerneltype == K_POLY, "Incorrect model kerneltype");
	mu_assert(model->gamma == 3.1, "Incorrect model gamma");
	mu_assert(model->coef == 2.1, "Incorrect model coef");
	mu_assert(model->degree == 1.1, "Incorrect model degree");
	// end test code //

	gensvm_free_model(model);
	gensvm_free_task(task);
	return NULL;
}

char *test_task_to_model_sigmoid()
{
	struct GenTask *task = gensvm_init_task();
	struct GenModel *model = gensvm_init_model();

	// start test code //
	task->weight_idx = 2;
	task->p = 1.3;
	task->kappa = 0.1;
	task->lambda = 1.4;
	task->epsilon = 5e-3;
	task->kerneltype = K_SIGMOID;
	task->gamma = 3.1;
	task->coef = 0.1;

	gensvm_task_to_model(task, model);

	mu_assert(model->weight_idx = 2, "Incorrect model weight_idx");
	mu_assert(model->p == 1.3, "Incorrect model p");
	mu_assert(model->kappa == 0.1, "Incorrect model kappa");
	mu_assert(model->lambda == 1.4, "Incorrect model lambda");
	mu_assert(model->epsilon == 5e-3, "Incorrect model epsilon");
	mu_assert(model->kerneltype == K_SIGMOID, "Incorrect model kerneltype");
	mu_assert(model->gamma == 3.1, "Incorrect model gamma");
	mu_assert(model->coef == 0.1, "Incorrect model coef");
	// end test code //

	gensvm_free_model(model);
	gensvm_free_task(task);
	return NULL;
}

char *test_copy_task_linear()
{
	struct GenTask *task = gensvm_init_task();
	struct GenTask *copy = NULL;
	struct GenData *train = gensvm_init_data();
	struct GenData *test = gensvm_init_data();

	// start test code //
	task->folds = 7;
	task->ID = 13;
	task->weight_idx = 2;
	task->p = 1.3;
	task->kappa = 0.1;
	task->lambda = 1.4;
	task->epsilon = 5e-3;
	task->kerneltype = K_LINEAR;
	task->train_data = train;
	task->test_data = test;
	task->performance = 11.11;

	copy = gensvm_copy_task(task);

	mu_assert(copy->folds == 7, "Incorrect copy folds");
	mu_assert(copy->ID == 13, "Incorrect copy ID");
	mu_assert(copy->weight_idx = 2, "Incorrect copy weight_idx");
	mu_assert(copy->p == 1.3, "Incorrect copy p");
	mu_assert(copy->kappa == 0.1, "Incorrect copy kappa");
	mu_assert(copy->lambda == 1.4, "Incorrect copy lambda");
	mu_assert(copy->epsilon == 5e-3, "Incorrect copy epsilon");
	mu_assert(copy->train_data == train, "Incorrect copy train data");
	mu_assert(copy->test_data == test, "Incorrect copy test data");
	mu_assert(copy->performance == 11.11, "Incorrect copy performance");
	mu_assert(copy->kerneltype == K_LINEAR, "Incorrect copy kerneltype");

	// end test code //
	gensvm_free_task(task);
	gensvm_free_task(copy);
	gensvm_free_data(train);
	gensvm_free_data(test);

	return NULL;
}

char *test_copy_task_rbf()
{
	struct GenTask *task = gensvm_init_task();
	struct GenTask *copy = NULL;
	struct GenData *train = gensvm_init_data();
	struct GenData *test = gensvm_init_data();

	// start test code //

	task->folds = 7;
	task->ID = 13;
	task->weight_idx = 2;
	task->p = 1.3;
	task->kappa = 0.1;
	task->lambda = 1.4;
	task->epsilon = 5e-3;
	task->train_data = train;
	task->test_data = test;
	task->performance = 11.11;
	task->kerneltype = K_RBF;
	task->gamma = 3.1;

	copy = gensvm_copy_task(task);

	mu_assert(copy->folds == 7, "Incorrect copy folds");
	mu_assert(copy->ID == 13, "Incorrect copy ID");
	mu_assert(copy->weight_idx = 2, "Incorrect copy weight_idx");
	mu_assert(copy->p == 1.3, "Incorrect copy p");
	mu_assert(copy->kappa == 0.1, "Incorrect copy kappa");
	mu_assert(copy->lambda == 1.4, "Incorrect copy lambda");
	mu_assert(copy->epsilon == 5e-3, "Incorrect copy epsilon");
	mu_assert(copy->train_data == train, "Incorrect copy train data");
	mu_assert(copy->test_data == test, "Incorrect copy test data");
	mu_assert(copy->performance == 11.11, "Incorrect copy performance");
	mu_assert(copy->kerneltype == K_RBF, "Incorrect copy kerneltype");
	mu_assert(copy->gamma == 3.1, "Incorrect copy gamma");

	// end test code //
	gensvm_free_task(copy);
	gensvm_free_task(task);
	gensvm_free_data(train);
	gensvm_free_data(test);

	return NULL;
}

char *test_copy_task_poly()
{
	struct GenTask *task = gensvm_init_task();
	struct GenTask *copy = NULL;
	struct GenData *train = gensvm_init_data();
	struct GenData *test = gensvm_init_data();

	// start test code //
	task->folds = 7;
	task->ID = 13;
	task->weight_idx = 2;
	task->p = 1.3;
	task->kappa = 0.1;
	task->lambda = 1.4;
	task->epsilon = 5e-3;
	task->train_data = train;
	task->test_data = test;
	task->performance = 11.11;
	task->kerneltype = K_POLY;
	task->gamma = 3.1;
	task->coef = 2.1;
	task->degree = 1.1;

	copy = gensvm_copy_task(task);

	mu_assert(copy->folds == 7, "Incorrect copy folds");
	mu_assert(copy->ID == 13, "Incorrect copy ID");
	mu_assert(copy->weight_idx = 2, "Incorrect copy weight_idx");
	mu_assert(copy->p == 1.3, "Incorrect copy p");
	mu_assert(copy->kappa == 0.1, "Incorrect copy kappa");
	mu_assert(copy->lambda == 1.4, "Incorrect copy lambda");
	mu_assert(copy->epsilon == 5e-3, "Incorrect copy epsilon");
	mu_assert(copy->train_data == train, "Incorrect copy train data");
	mu_assert(copy->test_data == test, "Incorrect copy test data");
	mu_assert(copy->performance == 11.11, "Incorrect copy performance");
	mu_assert(copy->kerneltype == K_POLY, "Incorrect copy kerneltype");
	mu_assert(copy->gamma == 3.1, "Incorrect copy gamma");
	mu_assert(copy->coef == 2.1, "Incorrect copy coef");
	mu_assert(copy->degree == 1.1, "Incorrect copy degree");

	// end test code //
	gensvm_free_task(task);
	gensvm_free_task(copy);
	gensvm_free_data(train);
	gensvm_free_data(test);

	return NULL;
}

char *test_copy_task_sigmoid()
{
	struct GenTask *task = gensvm_init_task();
	struct GenTask *copy = NULL;
	struct GenData *train = gensvm_init_data();
	struct GenData *test = gensvm_init_data();

	// start test code //
	task->folds = 7;
	task->ID = 13;
	task->weight_idx = 2;
	task->p = 1.3;
	task->kappa = 0.1;
	task->lambda = 1.4;
	task->epsilon = 5e-3;
	task->train_data = train;
	task->test_data = test;
	task->performance = 11.11;
	task->kerneltype = K_SIGMOID;
	task->gamma = 3.1;
	task->coef = 0.1;

	copy = gensvm_copy_task(task);

	mu_assert(copy->folds == 7, "Incorrect copy folds");
	mu_assert(copy->ID == 13, "Incorrect copy ID");
	mu_assert(copy->weight_idx = 2, "Incorrect copy weight_idx");
	mu_assert(copy->p == 1.3, "Incorrect copy p");
	mu_assert(copy->kappa == 0.1, "Incorrect copy kappa");
	mu_assert(copy->lambda == 1.4, "Incorrect copy lambda");
	mu_assert(copy->epsilon == 5e-3, "Incorrect copy epsilon");
	mu_assert(copy->train_data == train, "Incorrect copy train data");
	mu_assert(copy->test_data == test, "Incorrect copy test data");
	mu_assert(copy->performance == 11.11, "Incorrect copy performance");
	mu_assert(copy->kerneltype == K_SIGMOID, "Incorrect copy kerneltype");
	mu_assert(copy->gamma == 3.1, "Incorrect copy gamma");
	mu_assert(copy->coef == 0.1, "Incorrect copy coef");

	// end test code //
	gensvm_free_task(task);
	gensvm_free_task(copy);
	gensvm_free_data(train);
	gensvm_free_data(test);
	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_init_free_task);
	mu_run_test(test_task_to_model_linear);
	mu_run_test(test_task_to_model_rbf);
	mu_run_test(test_task_to_model_sigmoid);
	mu_run_test(test_task_to_model_poly);
	mu_run_test(test_copy_task_linear);
	mu_run_test(test_copy_task_rbf);
	mu_run_test(test_copy_task_sigmoid);
	mu_run_test(test_copy_task_poly);

	return NULL;
}

RUN_TESTS(all_tests);
