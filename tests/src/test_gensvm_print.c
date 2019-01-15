/**
 * @file test_gensvm_print.c
 * @author G.J.J. van den Burg
 * @date 2016-09-01
 * @brief Unit tests for gensvm_print.c functions
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
#include "gensvm_print.h"

extern FILE *GENSVM_OUTPUT_FILE;
extern FILE *GENSVM_ERROR_FILE;

char *test_note()
{
	FILE *fid = NULL;
	GENSVM_OUTPUT_FILE = fopen("./data/test_note_print.txt", "w");
	mu_assert(GENSVM_OUTPUT_FILE != NULL, "output file couldn't be opened");

	// start test code //
	note("This is some text.\n");
	note("This is %s with %.2f.\n", "formatted text", 1.231234);

	// close the output stream
	fclose(GENSVM_OUTPUT_FILE);
	GENSVM_OUTPUT_FILE = NULL;

	note("This is some more text.\n");
	note("This shouldn't appear in the output file.\n");

	char buffer[GENSVM_MAX_LINE_LENGTH];
	fid = fopen("./data/test_note_print.txt", "r");
	mu_assert(fid != NULL, "output file couldn't be opened for reading");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "This is some text.\n") == 0,
			"Line doesn't contain expected content (1)");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "This is formatted text with 1.23.\n") == 0,
			"Line doesn't contain expected content (2)");

	fclose(fid);

	// end test code //

	return NULL;
}

char *test_gensvm_error()
{
	FILE *fid = NULL;
	GENSVM_ERROR_FILE = fopen("./data/test_err_print.txt", "w");
	mu_assert(GENSVM_ERROR_FILE != NULL, "error file couldn't be opened");

	// start test code //
	gensvm_error("This is some text.\n");
	gensvm_error("This is %s with %.2f.\n", "formatted text", 1.231234);

	// close the output stream
	fclose(GENSVM_ERROR_FILE);
	GENSVM_ERROR_FILE = NULL;

	gensvm_error("This is some more text.\n");
	gensvm_error("This shouldn't appear in the output file.\n");

	char buffer[GENSVM_MAX_LINE_LENGTH];
	fid = fopen("./data/test_err_print.txt", "r");
	mu_assert(fid != NULL, "error file couldn't be opened for reading.");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	mu_assert(strcmp(buffer, "This is some text.\n") == 0,
			"Line doesn't contain expected content (1)");

	fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
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
	mu_run_test(test_gensvm_error);

	return NULL;
}

RUN_TESTS(all_tests);
