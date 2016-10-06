/**
 * @file test_gensvm_timer.c
 * @author Gertjan van den Burg
 * @date October, 2016
 * @brief Unit tests for gensvm_timer.c functions
 *
 */

#include "minunit.h"
#include "gensvm_timer.h"

char *test_timer()
{
	struct timespec start, stop;
	Timer(start);

	stop.tv_sec = start.tv_sec + 1;
	stop.tv_nsec = start.tv_nsec;

	mu_assert(gensvm_elapsed_time(&start, &stop) == 1.0,
			"Incorrect time computed");

	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_timer);

	return NULL;
}

RUN_TESTS(all_tests);
