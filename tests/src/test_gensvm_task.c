/**
 * @file test_gensvm_task.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Unit tests for gensvm_task.c functions
 */

#include "minunit.h"
#include "gensvm_task.h"

char *test_init_free_task()
{
	struct GenTask *task = gensvm_init_task();
	gensvm_free_task(task);
	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_init_free_task);

	return NULL;
}

RUN_TESTS(all_tests);
