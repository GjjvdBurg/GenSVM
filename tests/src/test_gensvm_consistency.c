/**
 * @file test_gensvm_consistency.c
 * @author G.J.J. van den Burg
 * @date 2016-10-24
 * @brief Unit tests for gensvm_consistency.c
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
#include "gensvm_consistency.h"

char *test_doublesort()
{
	double a = 1.0;
	double b = 2.0;

	mu_assert(doublesort(&b, &a) == true, "Incorrect doublesort (1)");
	mu_assert(doublesort(&a, &b) == false, "Incorrect doublesort (2)");
	mu_assert(doublesort(&a, &a) == false, "Incorrect doublesort (3)");

	return NULL;
}

char *test_percentile_1()
{
	double *values = Malloc(double, 1);
	values[0] = 0.1368311165400936;

	// start test code //
	mu_assert(fabs(gensvm_percentile(values, 1, 25.0) -
				0.1368311165400936) < 1e-16,
			"Incorrect percentile");
	// end test code //
	free(values);

	return NULL;
}

char *test_percentile()
{
	double *values = Malloc(double, 10);
	values[0] = 0.1368311165400936;
	values[1] = 0.0864373686918369;
	values[2] = 0.9959483430066688;
	values[3] = 0.2946638351338509;
	values[4] = 0.3535927892606028;
	values[5] = 0.5898175818278500;
	values[6] = 0.1769525979717794;
	values[7] = 0.3114487168265636;
	values[8] = 0.3895012665017124;
	values[9] = 0.3229492282960943;

	// start test code //
	mu_assert(fabs(gensvm_percentile(values, 10, 25.0) -
				0.176952597971779) < 1e-14,
			"Incorrect 25th percentile");
	mu_assert(fabs(gensvm_percentile(values, 10, 50.0) -
				0.317198972561329) < 1e-14,
			"Incorrect 50th percentile");
	mu_assert(fabs(gensvm_percentile(values, 10, 75.0) -
				0.389501266501712) < 1e-14,
			"Incorrect 75th percentile");
	mu_assert(fabs(gensvm_percentile(values, 10, 90.0) -
				0.792882962417259) < 1e-14,
			"Incorrect 90th percentile");
	// end test code //
	free(values);

	return NULL;
}

char *test_top_queue()
{
	int i, N = 10;
	struct GenQueue *q = gensvm_init_queue();
	q->tasks = Malloc(struct GenTask *, N);
	q->N = N;
	for (i=0; i<N; i++) {
		q->tasks[i] = gensvm_init_task();
		q->tasks[i]->ID = i+1;
	}

	q->tasks[0]->performance = 0.1368311165400936;
	q->tasks[1]->performance = 0.0864373686918369;
	q->tasks[2]->performance = 0.9959483430066688; //
	q->tasks[3]->performance = 0.2946638351338509;
	q->tasks[4]->performance = 0.3535927892606028;
	q->tasks[5]->performance = 0.5898175818278500; //
	q->tasks[6]->performance = 0.1769525979717794;
	q->tasks[7]->performance = 0.3114487168265636;
	q->tasks[8]->performance = 0.3895012665017124; //
	q->tasks[9]->performance = 0.3229492282960943;

	// start test code //

	// boundary should be determined at: 0.389501266501712
	struct GenQueue *nq = gensvm_top_queue(q, 75.0);
	mu_assert(nq->N == 3, "Incorrect size of top queue");

	// end test code //
	gensvm_free_queue(q);
	gensvm_free_queue(nq);

	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_doublesort);
	mu_run_test(test_percentile_1);
	mu_run_test(test_percentile);
	mu_run_test(test_top_queue);

	return NULL;
}

RUN_TESTS(all_tests);
