/**
 * @file test_gensvm_rand.c
 * @author G.J.J. van den Burg
 * @date 2018-04-17
 * @brief Unit tests for gensvm_rand.c functions
 *
 * @copyright
 Copyright 2018, G.J.J. van den Burg.

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
#include "gensvm_rand.h"

char *test_rand_1()
{
	srand(123);
	double a = rand();

	gensvm_srand(123);
	double b = gensvm_rand();

	mu_assert(a == b, "Random numbers unequal");

	return NULL;
}

char *test_rand_many()
{
	int i;

	srand(87431);
	gensvm_srand(87431);

	for (i=0; i<1000; i++)
		mu_assert(rand() == gensvm_rand(), "Random numbers unequal");

	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_rand_1);
	mu_run_test(test_rand_many);

	return NULL;
}

RUN_TESTS(all_tests);
