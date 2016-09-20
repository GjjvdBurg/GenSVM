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


char *test_get_tt_split()
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


char *all_tests()
{
	mu_suite_start();
	mu_run_test(test_make_cv_split_1);
	mu_run_test(test_make_cv_split_2);
	mu_run_test(test_get_tt_split);

	return NULL;
}

RUN_TESTS(all_tests);
