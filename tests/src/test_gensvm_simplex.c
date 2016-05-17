/**
 * @file test_gensvm_simplex.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Unit tests for gensvm_simplex.c functions
 */

#include "minunit.h"
#include "gensvm_simplex.h"

char *test_simplex_1()
{
	double *U = Calloc(double, 2*1);

	gensvm_simplex(2, U);

	mu_assert(matrix_get(U, 1, 0, 0) == -0.5, "U(0, 0) incorrect.");
	mu_assert(matrix_get(U, 1, 1, 0) == 0.5, "U(1, 0) incorrect.");

	free(U);

	return NULL;
}

char *test_simplex_2()
{
	double *U = Calloc(double, 4*3);

	gensvm_simplex(4, U);

	mu_assert(matrix_get(U, 3, 0, 0) == -0.5, "U(0, 0) incorrect.");
	mu_assert(matrix_get(U, 3, 1, 0) == 0.5, "U(1, 0) incorrect.");
	mu_assert(matrix_get(U, 3, 2, 0) == 0.0, "U(2, 0) incorrect.");
	mu_assert(matrix_get(U, 3, 3, 0) == 0.0, "U(3, 0) incorrect.");

	mu_assert(fabs(matrix_get(U, 3, 0, 1) - -0.5/sqrt(3)) < 1e-14,
		  	"U(0, 1) incorrect.");
	mu_assert(fabs(matrix_get(U, 3, 1, 1) - -0.5/sqrt(3)) < 1e-14,
		  	"U(1, 1) incorrect.");
	mu_assert(fabs(matrix_get(U, 3, 2, 1) - 1.0/sqrt(3)) < 1e-14,
		  	"U(2, 1) incorrect.");
	mu_assert(fabs(matrix_get(U, 3, 3, 1) - 0.0) < 1e-14,
		  	"U(3, 1) incorrect.");

	mu_assert(fabs(matrix_get(U, 3, 0, 2) - -1.0/sqrt(24)) < 1e-14,
		  	"U(0, 2) incorrect.");
	mu_assert(fabs(matrix_get(U, 3, 1, 2) - -1.0/sqrt(24)) < 1e-14,
		  	"U(1, 2) incorrect.");
	mu_assert(fabs(matrix_get(U, 3, 2, 2) - -1.0/sqrt(24)) < 1e-14,
		  	"U(2, 2) incorrect.");
	mu_assert(fabs(matrix_get(U, 3, 3, 2) - 3.0/sqrt(24)) < 1e-14,
		       	"U(3, 2) incorrect.");

	free(U);

	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_simplex_1);
	mu_run_test(test_simplex_2);

	return NULL;
}

RUN_TESTS(all_tests);
