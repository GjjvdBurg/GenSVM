/**
 * @file test_gensvm_debug.c
 * @author Gertjan van den Burg
 * @date September, 2016
 * @brief Unit tests for gensvm_io.c functions
 */

#include "minunit.h"
#include "gensvm_debug.h"

extern FILE *GENSVM_OUTPUT_FILE;

char *test_print_matrix()
{
	FILE *fid = NULL;
	GENSVM_OUTPUT_FILE = fopen("./data/test_debug_print.txt", "w");

	double *mat = Calloc(double, 3*2);
	matrix_set(mat, 2, 0, 0, -0.241053050258449);
	matrix_set(mat, 2, 0, 1, -0.599809408260836);
	matrix_set(mat, 2, 1, 0,  0.893318163305108);
	matrix_set(mat, 2, 1, 1, -0.344057630469285);
	matrix_set(mat, 2, 2, 0,  0.933948479216127);
	matrix_set(mat, 2, 2, 1, -0.474352026604967);

	// start test code //
	gensvm_print_matrix(mat, 3, 2);
	fclose(GENSVM_OUTPUT_FILE);

	char buffer[MAX_LINE_LENGTH];
	fid = fopen("./data/test_debug_print.txt", "r");

	fgets(buffer, MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "-0.241053 -0.599809\n") == 0,
		       	"Line doesn't contain expected content (0).\n");

	fgets(buffer, MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "+0.893318 -0.344058\n") == 0,
		       	"Line doesn't contain expected content (1).\n");

	fgets(buffer, MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "+0.933948 -0.474352\n") == 0,
		       	"Line doesn't contain expected content (2).\n");

	fclose(fid);
	// end test code //
	free(mat);

	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_print_matrix);

	return NULL;
}

RUN_TESTS(all_tests);
