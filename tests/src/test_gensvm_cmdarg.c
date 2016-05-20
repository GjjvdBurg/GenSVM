/**
 * @file test_gensvm_cmdarg.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Unit tests for gensvm_cmdarg.c functions
 */

#include "minunit.h"
#include "gensvm_cmdarg.h"

char *test_check_argv()
{
	char **argv = Malloc(char *, 3);
	argv[0] = "find";
	argv[1] = "this";
	argv[2] = "string";

	mu_assert(gensvm_check_argv(3, argv, "ring") == 2,
			"Failed on string 'ring'");
	mu_assert(gensvm_check_argv(3, argv, "hi") == 1,
			"Failed on string 'hi'");
	mu_assert(gensvm_check_argv(3, argv, "test") == 0,
			"Failed on string 'test'");
	free(argv);

	return NULL;
}

char *test_check_argv_eq()
{
	char **argv = Malloc(char *, 3);
	argv[0] = "find";
	argv[1] = "this";
	argv[2] = "string";

	mu_assert(gensvm_check_argv_eq(3, argv, "this") == 1,
			"Failed on string 'this'");
	mu_assert(gensvm_check_argv_eq(3, argv, "hi") == 0,
			"Failed on string 'hi'");
	mu_assert(gensvm_check_argv_eq(3, argv, "str") == 0,
			"Failed on string 'str'");
	free(argv);

	return NULL;
}


char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_check_argv);
	mu_run_test(test_check_argv_eq);

	return NULL;
}

RUN_TESTS(all_tests);
