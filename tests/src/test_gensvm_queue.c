/**
 * @file test_gensvm_queue.c
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Unit tests for gensvm_queue.c functions
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
#include "gensvm_queue.h"

char *test_init_free_queue()
{
	struct GenQueue *queue = gensvm_init_queue();
	gensvm_free_queue(queue);
	return NULL;
}

char *test_get_next_task()
{
	struct GenQueue *q = gensvm_init_queue();
	q->N = 3;
	q->i = 0;
	q->tasks = Calloc(struct GenTask *, 3);

	struct GenTask *task_1 = gensvm_init_task();
	struct GenTask *task_2 = gensvm_init_task();
	struct GenTask *task_3 = gensvm_init_task();
	task_1->ID = 1;
	task_2->ID = 2;
	task_3->ID = 3;

	q->tasks[0] = task_1;
	q->tasks[1] = task_2;
	q->tasks[2] = task_3;

	struct GenTask *t = NULL;
	t = get_next_task(q);
	mu_assert(t == task_1, "first task is not task_1");
	t = get_next_task(q);
	t = get_next_task(q);
	mu_assert(t == task_3, "third task is not task_3");
	t = get_next_task(q);
	mu_assert(t == NULL, "next task is not NULL");

	gensvm_free_queue(q);

	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_init_free_queue);
	mu_run_test(test_get_next_task);

	return NULL;
}

RUN_TESTS(all_tests);
