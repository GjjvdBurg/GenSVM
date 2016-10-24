/**
 * @file test_gensvm_cmdarg.c
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Unit tests for gensvm_cmdarg.c functions
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
