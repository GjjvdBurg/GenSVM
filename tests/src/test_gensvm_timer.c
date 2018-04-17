/**
 * @file test_gensvm_timer.c
 * @author G.J.J. van den Burg
 * @date 2016-10-01
 * @brief Unit tests for gensvm_timer.c functions
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
#include "gensvm_timer.h"

char *test_timer()
{
	GenTime start, stop;
	Timer(start);

	gensvm_sleep(1.0);

	Timer(stop);

	double elapsed = gensvm_elapsed_time(&start, &stop);
	mu_assert(elapsed > 0.99 && elapsed < 1.01, "Incorrect timing or sleeping.");

	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_timer);

	return NULL;
}

RUN_TESTS(all_tests);
