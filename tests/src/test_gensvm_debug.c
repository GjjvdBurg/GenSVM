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
	const char *filename = "./data/test_debug_dense.txt";
	GENSVM_OUTPUT_FILE = fopen(filename, "w");

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
	fid = fopen(filename, "r");

	fgets(buffer, MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "-0.241053 -0.599809\n") == 0,
		       	"Line doesn't contain expected content (0).");

	fgets(buffer, MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "+0.893318 -0.344058\n") == 0,
		       	"Line doesn't contain expected content (1).");

	fgets(buffer, MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "+0.933948 -0.474352\n") == 0,
		       	"Line doesn't contain expected content (2).");

	fclose(fid);
	// end test code //
	free(mat);

	return NULL;
}

char *test_print_sparse()
{
	FILE *fid = NULL;
	double *A = Calloc(double, 4*4);
	A[4] = 5;
	A[5] = 8;
	A[10] = 3;
	A[13] = 6;
	struct GenSparse *sp = gensvm_dense_to_sparse(A, 4, 4);
	const char *filename = "./data/test_debug_sparse.txt";
	GENSVM_OUTPUT_FILE = fopen(filename, "w");

	// start test code //
	gensvm_print_sparse(sp);
	fclose(GENSVM_OUTPUT_FILE);

	char buffer[MAX_LINE_LENGTH];
	fid = fopen(filename, "r");

	fgets(buffer, MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "Sparse Matrix:\n") == 0,
			"Line doesn't contain expected content (0).");

	fgets(buffer, MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "\tnnz = 4, rows = 4, cols = 4\n") == 0,
			"Line doesn't contain expected content (1).");

	fgets(buffer, MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "\tvalues = [ 5.000000, 8.000000, "
				"3.000000, 6.000000 ]\n") == 0,
			"Line doesn't contain expected content (2).");

	fgets(buffer, MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "\tIA = [ 0, 0, 2, 3, 4 ]\n") == 0,
			"Line doesn't contain expected content (3).");

	fgets(buffer, MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "\tJA = [ 0, 1, 2, 1 ]\n") == 0,
			"Line doesn't contain expected content (4).");

	fclose(fid);
	// end test code //
	gensvm_free_sparse(sp);
	free(A);

	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_print_matrix);
	mu_run_test(test_print_sparse);

	return NULL;
}

RUN_TESTS(all_tests);
