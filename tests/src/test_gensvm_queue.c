/**
 * @file test_gensvm_queue.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Unit tests for gensvm_queue.c functions
 */

#include "minunit.h"
#include "gensvm_task.h"
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
