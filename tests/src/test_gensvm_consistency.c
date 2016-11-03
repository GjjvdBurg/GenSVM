/**
 * @file test_gensvm_consistency.c
 * @author G.J.J. van den Burg
 * @date 2016-10-24
 * @brief Unit tests for gensvm_consistency.c
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

char *test_1()
{
	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_1);

	return NULL;
}

RUN_TESTS(all_tests);
