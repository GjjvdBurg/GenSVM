/**
 * @file test_gensvm_gridsearch.c
 * @author G.J.J. van den Burg
 * @date 2016-10-18
 * @brief Unit tests for gensvm_gridsearch.c
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
#include "gensvm_gridsearch.h"

char *test_fill_queue_nokernel()
{
	struct GenData *train_data = gensvm_init_data();
	struct GenData *test_data = gensvm_init_data();
	struct GenGrid *grid = gensvm_init_grid();

	grid->Np = 3;
	grid->Nl = 2;
	grid->Nk = 1;
	grid->Ne = 1;
	grid->Nw = 1;

	grid->ps = Calloc(double, grid->Np);
	grid->ps[0] = 1.0;
	grid->ps[1] = 1.5;
	grid->ps[2] = 2.0;

	grid->lambdas = Calloc(double, grid->Nl);
	grid->lambdas[0] = 1.0;
	grid->lambdas[1] = 5.0;

	grid->kappas = Calloc(double, grid->Nk);
	grid->kappas[0] = -0.99;

	grid->epsilons = Calloc(double, grid->Ne);
	grid->epsilons[0] = 1e-6;

	grid->weight_idxs = Calloc(double, grid->Nw);
	grid->weight_idxs[0] = 1;

	// start test code //
	struct GenQueue *q = gensvm_init_queue();
	gensvm_fill_queue(grid, q, train_data, test_data);

	mu_assert(q->N == 6, "Incorrect number of queue elements");

	int i;
	for (i=0; i<q->N; i++) {
		mu_assert(q->tasks[i]->ID == i, "Incorrect ID");
		mu_assert(q->tasks[i]->folds == 10, "Incorrect folds");
		mu_assert(q->tasks[i]->kerneltype == K_LINEAR, 
				"Incorrect kernel type");
		mu_assert(q->tasks[i]->train_data == train_data, 
				"Incorrect train_data");
		mu_assert(q->tasks[i]->test_data == test_data, 
				"Incorrect test data")
	}

	// test p value
	mu_assert(q->tasks[0]->p == 1.0, "Incorrect p at task 0");
	mu_assert(q->tasks[1]->p == 1.5, "Incorrect p at task 1");
	mu_assert(q->tasks[2]->p == 2.0, "Incorrect p at task 2");
	mu_assert(q->tasks[3]->p == 1.0, "Incorrect p at task 3");
	mu_assert(q->tasks[4]->p == 1.5, "Incorrect p at task 4");
	mu_assert(q->tasks[5]->p == 2.0, "Incorrect p at task 5");

	// test lambda value
	mu_assert(q->tasks[0]->lambda == 1.0, "Incorrect lambda at task 0");
	mu_assert(q->tasks[1]->lambda == 1.0, "Incorrect lambda at task 1");
	mu_assert(q->tasks[2]->lambda == 1.0, "Incorrect lambda at task 2");
	mu_assert(q->tasks[3]->lambda == 5.0, "Incorrect lambda at task 3");
	mu_assert(q->tasks[4]->lambda == 5.0, "Incorrect lambda at task 4");
	mu_assert(q->tasks[5]->lambda == 5.0, "Incorrect lambda at task 5");

	// test kappa value
	mu_assert(q->tasks[0]->kappa == -0.99, "Incorrect kappa at task 0");
	mu_assert(q->tasks[1]->kappa == -0.99, "Incorrect kappa at task 1");
	mu_assert(q->tasks[2]->kappa == -0.99, "Incorrect kappa at task 2");
	mu_assert(q->tasks[3]->kappa == -0.99, "Incorrect kappa at task 3");
	mu_assert(q->tasks[4]->kappa == -0.99, "Incorrect kappa at task 4");
	mu_assert(q->tasks[5]->kappa == -0.99, "Incorrect kappa at task 5");

	// test epsilon value
	mu_assert(q->tasks[0]->epsilon == 1e-6, "Incorrect epsilon at task 0");
	mu_assert(q->tasks[1]->epsilon == 1e-6, "Incorrect epsilon at task 1");
	mu_assert(q->tasks[2]->epsilon == 1e-6, "Incorrect epsilon at task 2");
	mu_assert(q->tasks[3]->epsilon == 1e-6, "Incorrect epsilon at task 3");
	mu_assert(q->tasks[4]->epsilon == 1e-6, "Incorrect epsilon at task 4");
	mu_assert(q->tasks[5]->epsilon == 1e-6, "Incorrect epsilon at task 5");

	gensvm_free_queue(q);
	// end test code //

	gensvm_free_data(train_data);
	gensvm_free_data(test_data);
	gensvm_free_grid(grid);

	return NULL;
}


char *test_kernel_changed()
{
	mu_test_missing();

	return NULL;
}

char *test_kernel_folds()
{
	mu_test_missing();

	return NULL;
}

char *test_train_queue()
{
	mu_test_missing();

	return NULL;
}

char *test_gridsearch_progress()
{
	mu_test_missing();

	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_fill_queue_nokernel);
	mu_run_test(test_kernel_changed);
	mu_run_test(test_kernel_folds);
	mu_run_test(test_train_queue);
	mu_run_test(test_gridsearch_progress);

	return NULL;
}

RUN_TESTS(all_tests);
