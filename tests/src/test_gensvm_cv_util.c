/**
 * @file test_gensvm_cv_util.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Unit tests for gensvm_cv_util.c functions
 */

#include "minunit.h"
#include "gensvm_cv_util.h"

char *test_make_cv_split_1()
{
	srand(0);
	int i, j;
	long N = 10;
	long folds = 4;
	long *cv_idx = Calloc(long, N);

	// start test code //
	gensvm_make_cv_split(N, folds, cv_idx);
	// check if the values are between [0, folds-1]
	for (i=0; i<N; i++)
		mu_assert(0 <= cv_idx[i] && cv_idx[i] < folds,
				"CV range incorrect.");

	// check there are N % folds big folds of size floor(N/folds) + 1
	// and the remaining are of size floor(N/folds)
	int sum;
	int is_big = 0,
	    is_small = 0;
	for (i=0; i<folds; i++) {
		sum = 0;
		for (j=0; j<N; j++) {
			if (cv_idx[j] == i) sum += 1;
		}
		if (sum == floor(N/folds) + 1)
			is_big++;
		else
			is_small++;
	}
	mu_assert(is_big == N % folds, "Incorrect number of big folds");
	mu_assert(is_small == folds - N % folds,
			"Incorrect number of small folds");

	// end test code //

	free(cv_idx);

	return NULL;
}

char *test_make_cv_split_2()
{
	srand(0);
	int i, j;
	long N = 101;
	long folds = 7;
	long *cv_idx = Calloc(long, N);

	// start test code //
	gensvm_make_cv_split(N, folds, cv_idx);
	// check if the values are between [0, folds-1]
	for (i=0; i<N; i++)
		mu_assert(0 <= cv_idx[i] && cv_idx[i] < folds,
				"CV range incorrect.");

	// check there are N % folds big folds of size floor(N/folds) + 1
	// and the remaining are of size floor(N/folds)
	int sum;
	int is_big = 0,
	    is_small = 0;
	for (i=0; i<folds; i++) {
		sum = 0;
		for (j=0; j<N; j++) {
			if (cv_idx[j] == i) sum += 1;
		}
		if (sum == floor(N/folds) + 1)
			is_big++;
		else
			is_small++;
	}
	mu_assert(is_big == N % folds, "Incorrect number of big folds");
	mu_assert(is_small == folds - N % folds,
			"Incorrect number of small folds");

	// end test code //

	free(cv_idx);

	return NULL;
}

char *test_get_tt_split_dense()
{
	struct GenData *full = gensvm_init_data();
	full->K = 3;
	full->n = 10;
	full->m = 2;
	full->r = 2;

	full->y = Calloc(long, full->n);
	full->y[0] = 1;
	full->y[1] = 2;
	full->y[2] = 3;
	full->y[3] = 1;
	full->y[4] = 2;
	full->y[5] = 3;
	full->y[6] = 1;
	full->y[7] = 2;
	full->y[8] = 3;
	full->y[9] = 1;

	full->RAW = Calloc(double, full->n * (full->m+1));
	matrix_set(full->RAW, full->m+1, 0, 1, 1.0);
	matrix_set(full->RAW, full->m+1, 0, 2, 1.0);
	matrix_set(full->RAW, full->m+1, 1, 1, 2.0);
	matrix_set(full->RAW, full->m+1, 1, 2, 2.0);
	matrix_set(full->RAW, full->m+1, 2, 1, 3.0);
	matrix_set(full->RAW, full->m+1, 2, 2, 3.0);
	matrix_set(full->RAW, full->m+1, 3, 1, 4.0);
	matrix_set(full->RAW, full->m+1, 3, 2, 4.0);
	matrix_set(full->RAW, full->m+1, 4, 1, 5.0);
	matrix_set(full->RAW, full->m+1, 4, 2, 5.0);
	matrix_set(full->RAW, full->m+1, 5, 1, 6.0);
	matrix_set(full->RAW, full->m+1, 5, 2, 6.0);
	matrix_set(full->RAW, full->m+1, 6, 1, 7.0);
	matrix_set(full->RAW, full->m+1, 6, 2, 7.0);
	matrix_set(full->RAW, full->m+1, 7, 1, 8.0);
	matrix_set(full->RAW, full->m+1, 7, 2, 8.0);
	matrix_set(full->RAW, full->m+1, 8, 1, 9.0);
	matrix_set(full->RAW, full->m+1, 8, 2, 9.0);
	matrix_set(full->RAW, full->m+1, 9, 1, 10.0);
	matrix_set(full->RAW, full->m+1, 9, 2, 10.0);
	full->Z = full->RAW;

	long *cv_idx = Calloc(long, full->n);
	cv_idx[0] = 1;
	cv_idx[1] = 0;
	cv_idx[2] = 1;
	cv_idx[3] = 0;
	cv_idx[4] = 1;
	cv_idx[5] = 2;
	cv_idx[6] = 3;
	cv_idx[7] = 2;
	cv_idx[8] = 3;
	cv_idx[9] = 2;

	struct GenData *train = gensvm_init_data();
	struct GenData *test = gensvm_init_data();

	// start test code //
	gensvm_get_tt_split(full, train, test, cv_idx, 0);

	mu_assert(train->n == 8, "train_n incorrect.");
	mu_assert(test->n == 2, "test_n incorrect.");

	mu_assert(train->m == 2, "train_m incorrect.");
	mu_assert(test->m == 2, "test_m incorrect.");

	mu_assert(train->K == 3, "train_K incorrect.");
	mu_assert(test->K == 3, "test_K incorrect.");

	mu_assert(train->y[0] == 1, "train y incorrect.");
	mu_assert(train->y[1] == 3, "train y incorrect.");
	mu_assert(train->y[2] == 2, "train y incorrect.");
	mu_assert(train->y[3] == 3, "train y incorrect.");
	mu_assert(train->y[4] == 1, "train y incorrect.");
	mu_assert(train->y[5] == 2, "train y incorrect.");
	mu_assert(train->y[6] == 3, "train y incorrect.");
	mu_assert(train->y[7] == 1, "train y incorrect.");

	mu_assert(test->y[0] == 2, "test y incorrect.");
	mu_assert(test->y[1] == 1, "test y incorrect.");

	mu_assert(matrix_get(train->RAW, train->m+1, 0, 0) == 0.0,
			"train RAW 0, 0 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 0, 1) == 1.0,
			"train RAW 0, 1 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 0, 2) == 1.0,
			"train RAW 0, 2 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 1, 0) == 0.0,
			"train RAW 1, 0 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 1, 1) == 3.0,
			"train RAW 1, 1 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 1, 2) == 3.0,
			"train RAW 1, 2 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 2, 0) == 0.0,
			"train RAW 2, 0 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 2, 1) == 5.0,
			"train RAW 2, 1 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 2, 2) == 5.0,
			"train RAW 2, 2 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 3, 0) == 0.0,
			"train RAW 3, 0 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 3, 1) == 6.0,
			"train RAW 3, 1 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 3, 2) == 6.0,
			"train RAW 3, 2 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 4, 0) == 0.0,
			"train RAW 4, 0 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 4, 1) == 7.0,
			"train RAW 4, 1 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 4, 2) == 7.0,
			"train RAW 4, 2 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 5, 0) == 0.0,
			"train RAW 5, 0 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 5, 1) == 8.0,
			"train RAW 5, 1 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 5, 2) == 8.0,
			"train RAW 5, 2 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 6, 0) == 0.0,
			"train RAW 6, 0 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 6, 1) == 9.0,
			"train RAW 6, 1 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 6, 2) == 9.0,
			"train RAW 6, 2 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 7, 0) == 0.0,
			"train RAW 7, 0 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 7, 1) == 10.0,
			"train RAW 7, 1 incorrect.");
	mu_assert(matrix_get(train->RAW, train->m+1, 7, 2) == 10.0,
			"train RAW 7, 2 incorrect.");

	mu_assert(matrix_get(test->RAW, train->m+1, 0, 0) == 0.0,
			"test RAW 0, 0 incorrect.");
	mu_assert(matrix_get(test->RAW, train->m+1, 0, 1) == 2.0,
			"test RAW 0, 1 incorrect.");
	mu_assert(matrix_get(test->RAW, train->m+1, 0, 2) == 2.0,
			"test RAW 0, 2 incorrect.");
	mu_assert(matrix_get(test->RAW, train->m+1, 1, 0) == 0.0,
			"test RAW 1, 0 incorrect.");
	mu_assert(matrix_get(test->RAW, train->m+1, 1, 1) == 4.0,
			"test RAW 1, 1 incorrect.");
	mu_assert(matrix_get(test->RAW, train->m+1, 1, 2) == 4.0,
			"test RAW 1, 2 incorrect.");

	// end test code //
	gensvm_free_data(full);
	gensvm_free_data(train);
	gensvm_free_data(test);
	free(cv_idx);

	return NULL;
}

char *test_get_tt_split_sparse()
{
	struct GenData *full = gensvm_init_data();
	full->K = 3;
	full->n = 10;
	full->m = 2;
	full->r = 2;

	full->y = Calloc(long, full->n);
	full->y[0] = 1;
	full->y[1] = 2;
	full->y[2] = 3;
	full->y[3] = 1;
	full->y[4] = 2;
	full->y[5] = 3;
	full->y[6] = 1;
	full->y[7] = 2;
	full->y[8] = 3;
	full->y[9] = 1;

	full->RAW = Calloc(double, full->n * (full->m+1));
	matrix_set(full->RAW, full->m+1, 0, 1, 1.0);
	matrix_set(full->RAW, full->m+1, 0, 2, 1.0);
	matrix_set(full->RAW, full->m+1, 1, 1, 2.0);
	matrix_set(full->RAW, full->m+1, 1, 2, 2.0);
	matrix_set(full->RAW, full->m+1, 2, 1, 3.0);
	matrix_set(full->RAW, full->m+1, 2, 2, 3.0);
	matrix_set(full->RAW, full->m+1, 3, 1, 4.0);
	matrix_set(full->RAW, full->m+1, 3, 2, 4.0);
	matrix_set(full->RAW, full->m+1, 4, 1, 5.0);
	matrix_set(full->RAW, full->m+1, 4, 2, 5.0);
	matrix_set(full->RAW, full->m+1, 5, 1, 6.0);
	matrix_set(full->RAW, full->m+1, 5, 2, 6.0);
	matrix_set(full->RAW, full->m+1, 6, 1, 7.0);
	matrix_set(full->RAW, full->m+1, 6, 2, 7.0);
	matrix_set(full->RAW, full->m+1, 7, 1, 8.0);
	matrix_set(full->RAW, full->m+1, 7, 2, 8.0);
	matrix_set(full->RAW, full->m+1, 8, 1, 9.0);
	matrix_set(full->RAW, full->m+1, 8, 2, 9.0);
	matrix_set(full->RAW, full->m+1, 9, 1, 10.0);
	matrix_set(full->RAW, full->m+1, 9, 2, 10.0);
	full->Z = full->RAW;

	// convert Z to a sparse matrix to test the sparse functions
	full->spZ = gensvm_dense_to_sparse(full->RAW, full->n, full->m+1);
	free(full->RAW);
	full->RAW = NULL;
	full->Z = NULL;

	long *cv_idx = Calloc(long, full->n);
	cv_idx[0] = 1;
	cv_idx[1] = 0;
	cv_idx[2] = 1;
	cv_idx[3] = 0;
	cv_idx[4] = 1;
	cv_idx[5] = 2;
	cv_idx[6] = 3;
	cv_idx[7] = 2;
	cv_idx[8] = 3;
	cv_idx[9] = 2;

	struct GenData *train = gensvm_init_data();
	struct GenData *test = gensvm_init_data();

	// start test code //
	gensvm_get_tt_split(full, train, test, cv_idx, 0);

	mu_assert(train->n == 8, "train_n incorrect.");
	mu_assert(test->n == 2, "test_n incorrect.");

	mu_assert(train->m == 2, "train_m incorrect.");
	mu_assert(test->m == 2, "test_m incorrect.");

	mu_assert(train->K == 3, "train_K incorrect.");
	mu_assert(test->K == 3, "test_K incorrect.");

	mu_assert(train->y[0] == 1, "train y incorrect.");
	mu_assert(train->y[1] == 3, "train y incorrect.");
	mu_assert(train->y[2] == 2, "train y incorrect.");
	mu_assert(train->y[3] == 3, "train y incorrect.");
	mu_assert(train->y[4] == 1, "train y incorrect.");
	mu_assert(train->y[5] == 2, "train y incorrect.");
	mu_assert(train->y[6] == 3, "train y incorrect.");
	mu_assert(train->y[7] == 1, "train y incorrect.");

	mu_assert(test->y[0] == 2, "test y incorrect.");
	mu_assert(test->y[1] == 1, "test y incorrect.");

	// check the train GenSparse struct
	mu_assert(train->spZ->nnz == 16, "train nnz incorrect");
	mu_assert(train->spZ->n_row == 8, "train n_row incorrect");
	mu_assert(train->spZ->n_col == 3, "train n_col incorrect");

	mu_assert(train->spZ->values[0] == 1.0, "Wrong train value at 0");
	mu_assert(train->spZ->values[1] == 1.0, "Wrong train value at 1");
	mu_assert(train->spZ->values[2] == 3.0, "Wrong train value at 2");
	mu_assert(train->spZ->values[3] == 3.0, "Wrong train value at 3");
	mu_assert(train->spZ->values[4] == 5.0, "Wrong train value at 4");
	mu_assert(train->spZ->values[5] == 5.0, "Wrong train value at 5");
	mu_assert(train->spZ->values[6] == 6.0, "Wrong train value at 6");
	mu_assert(train->spZ->values[7] == 6.0, "Wrong train value at 7");
	mu_assert(train->spZ->values[8] == 7.0, "Wrong train value at 8");
	mu_assert(train->spZ->values[9] == 7.0, "Wrong train value at 9");
	mu_assert(train->spZ->values[10] == 8.0, "Wrong train value at 10");
	mu_assert(train->spZ->values[11] == 8.0, "Wrong train value at 11");
	mu_assert(train->spZ->values[12] == 9.0, "Wrong train value at 12");
	mu_assert(train->spZ->values[13] == 9.0, "Wrong train value at 13");
	mu_assert(train->spZ->values[14] == 10.0, "Wrong train value at 14");
	mu_assert(train->spZ->values[15] == 10.0, "Wrong train value at 15");

	mu_assert(train->spZ->ia[0] == 0, "Wrong train ia at 0");
	mu_assert(train->spZ->ia[1] == 2, "Wrong train ia at 1");
	mu_assert(train->spZ->ia[2] == 4, "Wrong train ia at 2");
	mu_assert(train->spZ->ia[3] == 6, "Wrong train ia at 3");
	mu_assert(train->spZ->ia[4] == 8, "Wrong train ia at 4");
	mu_assert(train->spZ->ia[5] == 10, "Wrong train ia at 5");
	mu_assert(train->spZ->ia[6] == 12, "Wrong train ia at 6");
	mu_assert(train->spZ->ia[7] == 14, "Wrong train ia at 7");
	mu_assert(train->spZ->ia[8] == 16, "Wrong train ia at 8");

	mu_assert(train->spZ->ja[0] == 1, "Wrong train ja at 0");
	mu_assert(train->spZ->ja[1] == 2, "Wrong train ja at 1");
	mu_assert(train->spZ->ja[2] == 1, "Wrong train ja at 2");
	mu_assert(train->spZ->ja[3] == 2, "Wrong train ja at 3");
	mu_assert(train->spZ->ja[4] == 1, "Wrong train ja at 4");
	mu_assert(train->spZ->ja[5] == 2, "Wrong train ja at 5");
	mu_assert(train->spZ->ja[6] == 1, "Wrong train ja at 6");
	mu_assert(train->spZ->ja[7] == 2, "Wrong train ja at 7");
	mu_assert(train->spZ->ja[8] == 1, "Wrong train ja at 8");
	mu_assert(train->spZ->ja[9] == 2, "Wrong train ja at 9");
	mu_assert(train->spZ->ja[10] == 1, "Wrong train ja at 10");
	mu_assert(train->spZ->ja[11] == 2, "Wrong train ja at 11");
	mu_assert(train->spZ->ja[12] == 1, "Wrong train ja at 12");
	mu_assert(train->spZ->ja[13] == 2, "Wrong train ja at 13");
	mu_assert(train->spZ->ja[14] == 1, "Wrong train ja at 14");
	mu_assert(train->spZ->ja[15] == 2, "Wrong train ja at 15");

	// check the test GenSparse struct
	mu_assert(test->spZ->nnz == 4, "test nnz incorrect");
	mu_assert(test->spZ->n_row == 2, "test n_row incorrect");
	mu_assert(test->spZ->n_col == 3, "test n_col incorrect");

	mu_assert(test->spZ->values[0] == 2.0, "Wrong test value at 0");
	mu_assert(test->spZ->values[1] == 2.0, "Wrong test value at 1");
	mu_assert(test->spZ->values[2] == 4.0, "Wrong test value at 2");
	mu_assert(test->spZ->values[3] == 4.0, "Wrong test value at 3");

	mu_assert(test->spZ->ia[0] == 0, "Wrong test ia at 0");
	mu_assert(test->spZ->ia[1] == 2, "Wrong test ia at 1");
	mu_assert(test->spZ->ia[2] == 4, "Wrong test ia at 2");

	mu_assert(test->spZ->ja[0] == 1, "Wrong test ja at 0");
	mu_assert(test->spZ->ja[1] == 2, "Wrong test ja at 1");
	mu_assert(test->spZ->ja[2] == 1, "Wrong test ja at 2");
	mu_assert(test->spZ->ja[3] == 2, "Wrong test ja at 3");

	// end test code //
	gensvm_free_data(full);
	gensvm_free_data(train);
	gensvm_free_data(test);
	free(cv_idx);

	return NULL;
}

char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_make_cv_split_1);
	mu_run_test(test_make_cv_split_2);
	mu_run_test(test_get_tt_split_dense);
	mu_run_test(test_get_tt_split_sparse);

	return NULL;
}

RUN_TESTS(all_tests);
