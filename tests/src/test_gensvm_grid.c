/**
 * @file test_gensvm_grid.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Unit tests for gensvm_grid.c functions
 */

#include "minunit.h"
#include "gensvm_grid.h"

char *test_init_free_grid()
{
	struct GenGrid *grid = gensvm_init_grid();
	gensvm_free_grid(grid);
	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_init_free_grid);

	return NULL;
}

RUN_TESTS(all_tests);
