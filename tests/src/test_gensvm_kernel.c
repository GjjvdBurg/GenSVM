/**
 *@file test_gensvm_kernel.c
 *@author G.J.J. van den Burg
 *@date 2016-11-09
 *@brief Unit tests for gensvm_kernel.c
 *
 *@copyright
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
#include "gensvm_kernel.h"

char *test_dot_rbf()
{
	double dot;
	double *a = Malloc(double, 5);
	double *b = Malloc(double, 5);
	double *kernelparam = Malloc(double, 1);

	a[0] = 0.5203363837176203;
	a[1] = 0.3860628599460129;
	a[2] = 0.3592536954640216;
	a[3] = 0.6824659760765744;
	a[4] = 0.5390520090020700;

	b[0] = 0.1782643262351465;
	b[1] = 0.0314270210724957;
	b[2] = 0.5887219369641497;
	b[3] = 0.7710042954911620;
	b[4] = 0.8805451245738238;

	// start test code //
	kernelparam[0] = 1.0;
	dot = gensvm_kernel_dot_rbf(a, b, kernelparam, 5);
	mu_assert(fabs(dot - 0.657117701533133) < 1e-14, "Incorrect dot (1)");

	kernelparam[0] = 5.0;
	dot = gensvm_kernel_dot_rbf(a, b, kernelparam, 5);
	mu_assert(fabs(dot - 0.122522495044048) < 1e-14, "Incorrect dot (2)");

	// end test code //
	free(a);
	free(b);
	free(kernelparam);

	return NULL;
}

char *test_dot_poly()
{
	double dot;
	double *a = Malloc(double, 5);
	double *b = Malloc(double, 5);
	double *kernelparam = Malloc(double, 3);

	a[0] = 0.5203363837176203;
	a[1] = 0.3860628599460129;
	a[2] = 0.3592536954640216;
	a[3] = 0.6824659760765744;
	a[4] = 0.5390520090020700;

	b[0] = 0.1782643262351465;
	b[1] = 0.0314270210724957;
	b[2] = 0.5887219369641497;
	b[3] = 0.7710042954911620;
	b[4] = 0.8805451245738238;

	// start test code //
	kernelparam[0] = 1.0;
	kernelparam[1] = 1.0;
	kernelparam[2] = 1.0;
	dot = gensvm_kernel_dot_poly(a, b, kernelparam, 5);
	mu_assert(fabs(dot - 2.31723456944910) < 1e-14, "Incorrect dot (1)");

	kernelparam[0] = 1.5;
	kernelparam[1] = 2.5;
	kernelparam[2] = 3.5;
	dot = gensvm_kernel_dot_poly(a, b, kernelparam, 5);
	mu_assert(fabs(dot - 189.6989652572890179) < 1e-14, "Incorrect dot (2)");

	// end test code //
	free(a);
	free(b);
	free(kernelparam);
	return NULL;
}

char *test_dot_sigmoid()
{
	double dot;
	double *a = Malloc(double, 5);
	double *b = Malloc(double, 5);
	double *kernelparam = Malloc(double, 2);

	a[0] = 0.5203363837176203;
	a[1] = 0.3860628599460129;
	a[2] = 0.3592536954640216;
	a[3] = 0.6824659760765744;
	a[4] = 0.5390520090020700;

	b[0] = 0.1782643262351465;
	b[1] = 0.0314270210724957;
	b[2] = 0.5887219369641497;
	b[3] = 0.7710042954911620;
	b[4] = 0.8805451245738238;

	// start test code //
	kernelparam[0] = 1.0;
	kernelparam[1] = 1.0;
	dot = gensvm_kernel_dot_sigmoid(a, b, kernelparam, 5);
	mu_assert(fabs(dot - 0.9807642810850747) < 1e-14, "Incorrect dot (1)");

	kernelparam[0] = 1.5;
	kernelparam[1] = 2.5;
	dot = gensvm_kernel_dot_sigmoid(a, b, kernelparam, 5);
	mu_assert(fabs(dot - 0.9997410009167159) < 1e-14, "Incorrect dot (2)");

	// end test code //
	free(a);
	free(b);
	free(kernelparam);

	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_dot_rbf);
	mu_run_test(test_dot_poly);
	mu_run_test(test_dot_sigmoid);

	return NULL;
}

RUN_TESTS(all_tests);
