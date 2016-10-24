/**
 * @file test_gensvm_sparse.c
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Unit tests for gensvm_sparse.c functions
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
#include "gensvm_sparse.h"

char *test_init_free_sparse()
{
	struct GenSparse *sp = gensvm_init_sparse();

	mu_assert(sp->nnz == 0, "nnz not initialized correctly");
	mu_assert(sp->n_row == 0, "n not initialized correctly");
	mu_assert(sp->n_col == 0, "m not initialized correctly");
	mu_assert(sp->values == NULL, "values not initialized correctly");
	mu_assert(sp->ia == NULL, "ia not initialized correctly");
	mu_assert(sp->ja == NULL, "ja not initialized correctly");

	gensvm_free_sparse(sp);

	return NULL;
}

char *test_count_nnz()
{
	double *A = Calloc(double, 3*2);
	A[0] = 1.0;
	A[3] = 1.0;
	mu_assert(gensvm_count_nnz(A, 3, 2) == 2, "Incorrect nnz (1)");

	A[1] = 3.0;
	A[4] = 1e-20;
	mu_assert(gensvm_count_nnz(A, 3, 2) == 4, "Incorrect nnz (2)");

	free(A);

	return NULL;
}

char *test_gensvm_could_sparse()
{
	double *A = Calloc(double, 5*2);
	A[0] = 1.0;

	mu_assert(gensvm_could_sparse(A, 5, 2) == true, 
			"Incorrect could sparse (1)");
	A[1] = -1.0;
	mu_assert(gensvm_could_sparse(A, 5, 2) == false, 
			"Incorrect could sparse (2)");

	free(A);
	return NULL;
}

char *test_dense_to_sparse()
{
	double *A = Calloc(double, 4*4);
	A[4] = 5;
	A[5] = 8;
	A[10] = 3;
	A[13] = 6;

	struct GenSparse *sp = gensvm_dense_to_sparse(A, 4, 4);
	mu_assert(sp->nnz == 4, "Incorrect nnz");
	mu_assert(sp->n_row == 4, "Incorrect n_row");
	mu_assert(sp->n_col == 4, "Incorrect n_col");

	mu_assert(sp->values[0] == 5.0, "Incorrect value at 0");
	mu_assert(sp->values[1] == 8.0, "Incorrect value at 1");
	mu_assert(sp->values[2] == 3.0, "Incorrect value at 2");
	mu_assert(sp->values[3] == 6.0, "Incorrect value at 3");

	mu_assert(sp->ia[0] == 0, "Incorrect ia at 0");
	mu_assert(sp->ia[1] == 0, "Incorrect ia at 1");
	mu_assert(sp->ia[2] == 2, "Incorrect ia at 2");
	mu_assert(sp->ia[3] == 3, "Incorrect ia at 3");
	mu_assert(sp->ia[4] == 4, "Incorrect ia at 4");

	mu_assert(sp->ja[0] == 0, "Incorrect ja at 0");
	mu_assert(sp->ja[1] == 1, "Incorrect ja at 1");
	mu_assert(sp->ja[2] == 2, "Incorrect ja at 2");
	mu_assert(sp->ja[3] == 1, "Incorrect ja at 3");

	gensvm_free_sparse(sp);
	free(A);

	return NULL;
}

char *test_sparse_to_dense()
{
	double *A = Calloc(double, 4*4);
	A[4] = 5;
	A[5] = 8;
	A[10] = 3;
	A[13] = 6;

	struct GenSparse *sp = gensvm_dense_to_sparse(A, 4, 4);

	double *B = gensvm_sparse_to_dense(sp);
	int i;
	for (i=0; i<4*4; i++)
		mu_assert(B[i] == A[i], "Incorrect element of B");

	gensvm_free_sparse(sp);
	free(A);
	free(B);

	return NULL;
}

char *all_tests()
{
	mu_suite_start();

	mu_run_test(test_init_free_sparse);
	mu_run_test(test_count_nnz);
	mu_run_test(test_gensvm_could_sparse);
	mu_run_test(test_dense_to_sparse);
	mu_run_test(test_sparse_to_dense);

	return NULL;
}

RUN_TESTS(all_tests);
