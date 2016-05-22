/**
 * @file test_gensvm_strutil.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Unit tests for gensvm_strutil.c functions
 */

#include "minunit.h"
#include "gensvm_strutil.h"

char *test_str_startswith()
{
	mu_assert(str_startswith("test this string", "test") == true,
			"startswith first test failed.");
	mu_assert(str_startswith("test this string", "word") == false,
			"startswith second test failed.");
	return NULL;
}

char *test_str_endswith()
{
	mu_assert(str_endswith("this is a string", "string") == true,
			"endswith first test failed.");
	mu_assert(str_endswith("this is a string", "word") == false,
			"endswith second test failed.");
	return NULL;
}

char *test_get_line()
{
	char *fname = "test_get_line.txt";
	// prepare the test file
	FILE *fid = fopen(fname, "w");
	if (fid == NULL) return "couldn't open test file for get_line";
	fprintf(fid, "first line\n");
	fprintf(fid, "second line\n");
	fclose(fid);
	// end preparation

	// start of test code
	char buffer[MAX_LINE_LENGTH];
	char *retval = NULL;
	fid = fopen(fname, "r");
	retval = get_line(fid, fname, buffer);
	mu_assert(retval == buffer, "return value not buffer (1)");
	retval = get_line(fid, fname, buffer);
	mu_assert(retval == buffer, "return value not buffer (2)");
	retval = get_line(fid, fname, buffer);
	mu_assert(retval == NULL, "return value not NULL");

	fclose(fid);
	// end of test code

	// cleanup
	remove(fname);

	return NULL;
}

char *test_next_line()
{
	char *fname = "test_next_line.txt";
	// prepare the test file
	FILE *fid = fopen(fname, "w");
	if (fid == NULL) return "couldn't open test file for next_line";
	fprintf(fid, "first line\n");
	fprintf(fid, "second line\n");
	fprintf(fid, "third line\n");
	fprintf(fid, "fourth line\n");
	fclose(fid);
	// end preparation

	// start of test code
	char line[60];
	fid = fopen(fname, "r");
	next_line(fid, fname);
	fgets(line, 60, fid);
	mu_assert(strcmp(line, "second line\n") == 0, "second line unequal");
	next_line(fid, fname);
	fgets(line, 60, fid);
	mu_assert(strcmp(line, "fourth line\n") == 0, "fourth line unequal");
	fclose(fid);
	// end of test code

	// cleanup
	remove(fname);

	return NULL;
}

char *test_get_fmt_double()
{
	char *fname = "test_get_fmt_double.txt";
	// prepare the test file
	FILE *fid = fopen(fname, "w");
	if (fid == NULL) return "couldn't open test file for get_fmt_double";
	fprintf(fid, "double = 3.1416\n");
	fprintf(fid, "another double = 42.42\n");
	fprintf(fid, "line without a double\n");
	fclose(fid);

	// start test code //
	double value = -1.0;
	fid = fopen(fname, "r");

	value = get_fmt_double(fid, fname, "double = %lf");
	mu_assert(value == 3.1416, "first value wrong");

	value = get_fmt_double(fid, fname, "another double = %lf");
	mu_assert(value == 42.42, "second value wrong");

	value = get_fmt_double(fid, fname, "line without a double");
	mu_assert(isnan(value), "third value wrong");

	fclose(fid);
	// end test code //

	// cleanup
	remove(fname);

	return NULL;
}

char *test_get_fmt_long()
{
	char *fname = "test_get_fmt_long.txt";
	// prepare the test file
	FILE *fid = fopen(fname, "w");
	if (fid == NULL) return "couldn't open test file for get_fmt_double";
	fprintf(fid, "long = 53161\n");
	fprintf(fid, "another long = 1212\n");
	fprintf(fid, "line without a long\n");
	fclose(fid);

	// start test code //
	long value = -1;
	fid = fopen(fname, "r");

	value = get_fmt_long(fid, fname, "long = %li");
	mu_assert(value == 53161, "first value wrong");

	value = get_fmt_long(fid, fname, "another long = %li");
	mu_assert(value == 1212, "second value wrong");

	value = get_fmt_long(fid, fname, "line without a long");
	// unfortunately we can't test this properly, a warning will be 
	// generated but the return value will be 0, which can be a valid 
	// value. In general therefore, the real test is to see if a warning 
	// is printed to stderr.
	mu_assert(value == 0, "third value wrong")

		fclose(fid);
	// end test code //

	// cleanup
	remove(fname);

	return NULL;
}

char *test_all_doubles_str()
{
	char *fname = "test_all_doubles_str.txt";
	// prepare the test file
	FILE *fid = fopen(fname, "w");
	if (fid == NULL) return "couldn't open test file for all_doubles_str";
	fprintf(fid, "1.0 2.0 3.0 4.0\n");
	fprintf(fid, "9.1 8.2 7.3\n");
	fprintf(fid, "offset 1.0 2.0\n");
	fclose(fid);

	// start test code //
	char buffer[60];
	double *values = Calloc(double, 10);
	fid = fopen(fname, "r");
	int nr = -1;

	fgets(buffer, 60, fid);
	nr = all_doubles_str(buffer, 0, values);
	mu_assert(nr == 4, "incorrect number of doubles (1)");
	mu_assert(values[0] == 1.0, "incorrect first value (1)");
	mu_assert(values[1] == 2.0, "incorrect second value (1)");
	mu_assert(values[2] == 3.0, "incorrect third value (1)");
	mu_assert(values[3] == 4.0, "incorrect fourth value (1)");

	fgets(buffer, 60, fid);
	nr = all_doubles_str(buffer, 0, values);
	mu_assert(nr == 3, "incorrect number of doubles (2)");
	mu_assert(values[0] == 9.1, "incorrect first value (2)");
	mu_assert(values[1] == 8.2, "incorrect second value (2)");
	mu_assert(values[2] == 7.3, "incorrect third value (2)");

	fgets(buffer, 60, fid);
	nr = all_doubles_str(buffer, 7, values);
	mu_assert(nr == 2, "incorrect number of doubles (3)");
	mu_assert(values[0] == 1.0, "incorrect first value (3)");
	mu_assert(values[1] == 2.0, "incorrect second value (3)");
	// end test code //

	// cleanup
	remove(fname);
	free(values);

	return NULL;
}

char *test_all_longs_str()
{
	char *fname = "test_all_longs_str.txt";
	// prepare the test file
	FILE *fid = fopen(fname, "w");
	if (fid == NULL) return "couldn't open test file for all_longs_str";
	fprintf(fid, "1 2 3 4\n");
	fprintf(fid, "91 82 73\n");
	fprintf(fid, "offset 10 20\n");
	fclose(fid);

	// start test code //
	char buffer[60];
	long *values = Calloc(long, 10);
	fid = fopen(fname, "r");
	int nr = -1;

	fgets(buffer, 60, fid);
	nr = all_longs_str(buffer, 0, values);
	mu_assert(nr == 4, "incorrect number of longs (1)");
	mu_assert(values[0] == 1, "incorrect first value (1)");
	mu_assert(values[1] == 2, "incorrect second value (1)");
	mu_assert(values[2] == 3, "incorrect third value (1)");
	mu_assert(values[3] == 4, "incorrect fourth value (1)");

	fgets(buffer, 60, fid);
	nr = all_longs_str(buffer, 0, values);
	mu_assert(nr == 3, "incorrect number of longs (2)");
	mu_assert(values[0] == 91, "incorrect first value (2)");
	mu_assert(values[1] == 82, "incorrect second value (2)");
	mu_assert(values[2] == 73, "incorrect third value (2)");

	fgets(buffer, 60, fid);
	nr = all_longs_str(buffer, 7, values);
	mu_assert(nr == 2, "incorrect number of longs (3)");
	mu_assert(values[0] == 10, "incorrect first value (3)");
	mu_assert(values[1] == 20, "incorrect second value (3)");
	// end test code //

	// cleanup
	remove(fname);
	free(values);

	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_str_startswith);
	mu_run_test(test_str_endswith);
	mu_run_test(test_get_line);
	mu_run_test(test_next_line);
	mu_run_test(test_get_fmt_double);
	mu_run_test(test_get_fmt_long);
	mu_run_test(test_all_doubles_str);
	mu_run_test(test_all_longs_str);

	return NULL;
}

RUN_TESTS(all_tests);
