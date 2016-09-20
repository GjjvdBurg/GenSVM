/**
 * @file test_gensvm_print.c
 * @author Gertjan van den Burg
 * @date September, 2016
 * @brief Unit tests for gensvm_print.c functions
 */

#include "minunit.h"
#include "gensvm_print.h"

extern FILE *GENSVM_OUTPUT_FILE;
extern FILE *GENSVM_ERROR_FILE;

char *test_note()
{
	FILE *fid = NULL;
	GENSVM_OUTPUT_FILE = fopen("./data/test_note_print.txt", "w");

	// start test code //
	note("This is some text.\n");
	note("This is %s with %.2f.\n", "formatted text", 1.231234);

	// close the output stream
	fclose(GENSVM_OUTPUT_FILE);
	GENSVM_OUTPUT_FILE = NULL;

	note("This is some more text.\n");
	note("This shouldn't appear in the output file.\n");

	char buffer[MAX_LINE_LENGTH];
	fid = fopen("./data/test_note_print.txt", "r");

	fgets(buffer, MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "This is some text.\n") == 0,
			"Line doesn't contain expected content (1)");

	fgets(buffer, MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "This is formatted text with 1.23.\n") == 0,
			"Line doesn't contain expected content (2)");

	fclose(fid);

	// end test code //

	return NULL;
}

char *test_err()
{
	FILE *fid = NULL;
	GENSVM_ERROR_FILE = fopen("./data/test_err_print.txt", "w");

	// start test code //
	err("This is some text.\n");
	err("This is %s with %.2f.\n", "formatted text", 1.231234);

	// close the output stream
	fclose(GENSVM_ERROR_FILE);
	GENSVM_ERROR_FILE = NULL;

	err("This is some more text.\n");
	err("This shouldn't appear in the output file.\n");

	char buffer[MAX_LINE_LENGTH];
	fid = fopen("./data/test_err_print.txt", "r");

	fgets(buffer, MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "This is some text.\n") == 0,
			"Line doesn't contain expected content (1)");

	fgets(buffer, MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "This is formatted text with 1.23.\n") == 0,
			"Line doesn't contain expected content (2)");

	fclose(fid);

	// end test code //
	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_note);
	mu_run_test(test_err);

	return NULL;
}

RUN_TESTS(all_tests);
