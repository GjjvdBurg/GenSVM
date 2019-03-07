/**
 * @file gensvm_cv_util.c
 * @author G.J.J. van den Burg
 * @date 2014-01-07
 * @brief Functions for cross validation
 *
 * @details
 * This file contains functions for performing cross validation. The funtion
 * gensvm_make_cv_split() creates a cross validation vector for non-stratified
 * cross validation. The function gensvm_get_tt_split() creates a train and
 * test dataset from a given dataset and a pre-determined CV partition vector.
 * See individual function documentation for details.
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

#include "gensvm_cv_util.h"

/**
 * @brief Create a cross validation split vector
 *
 * @details
 * A pre-allocated vector of length N is created which can be used to define
 * cross validation splits. The folds are contain between
 * @f$ \lfloor N / folds \rfloor @f$ and @f$ \lceil N / folds \rceil @f$
 * instances. An instance is mapped to a partition randomly until all folds
 * contain @f$ N \% folds @f$ instances. The zero fold then contains
 * @f$ N / folds + N \% folds @f$ instances. These remaining @f$ N \% folds @f$
 * instances are then distributed over the first @f$ N \% folds @f$ folds.
 *
 * @param[in] 		N 	number of instances
 * @param[in] 		folds 	number of folds
 * @param[in,out] 	cv_idx 	array of size N which contains the fold index
 * 				for each observation on exit
 *
 */
void gensvm_make_cv_split(long N, long folds, long *cv_idx)
{
	long i, j, idx;

	for (i=0; i<N; i++)
		cv_idx[i] = 0;

	long big_folds = N%folds;
	long small_fold_size = N/folds;

	j = 0;
	for (i=0; i<small_fold_size*folds; i++)
		while (1) {
			idx = gensvm_rand() % N;
			if (cv_idx[idx] == 0) {
				cv_idx[idx] = j;
				j++;
				j%=folds;
				break;
			}
		}
	j = 0;
	i = 0;
	while (i < big_folds) {
		if (cv_idx[j] == 0) {
			cv_idx[j] = i++;
		}
		j++;
	}
}

/**
 * @brief Wrapper around sparse/dense versions of this function
 *
 * @details
 * This function tests if the data in the full structure is stored in a
 * dense matrix format or not, and calls gensvm_get_tt_split_dense() or
 * gensvm_get_tt_split_sparse() accordingly.
 *
 * @sa
 * gensvm_get_tt_split_dense(), gensvm_get_tt_split_sparse()
 *
 * @param[in] 		full 		a GenData structure for the entire
 * 					dataset
 * @param[in,out] 	train 		an initialized GenData structure which
 * 					on exit contains the training dataset
 * @param[in,out] 	test 		an initialized GenData structure which
 * 					on exit contains the test dataset
 * @param[in] 		cv_idx 		a vector of cv partitions created by
 * 					gensvm_make_cv_split()
 * @param[in] 		fold_idx 	index of the fold which becomes the
 * 					test dataset
 */
void gensvm_get_tt_split(struct GenData *full,
		struct GenData *train, struct GenData *test,
		long *cv_idx, long fold_idx)
{
	if (full->Z == NULL)
		gensvm_get_tt_split_sparse(full, train, test,
				cv_idx, fold_idx);
	else
		gensvm_get_tt_split_dense(full, train, test,
				cv_idx, fold_idx);
}

/**
 * @brief Create train and test datasets for a CV split with dense data
 *
 * @details
 * Given a GenData structure for the full dataset, a previously created
 * cross validation split vector and a fold index, a training and test dataset
 * are created. It is assumed here that the data is stored as a dense matrix,
 * and that the train and test data should also be stored as a dense matrix.
 *
 * @sa
 * gensvm_get_tt_split_sparse(), gensvm_get_tt_split()
 *
 * @param[in] 		full 		a GenData structure for the entire
 * 					dataset
 * @param[in,out] 	train 		an initialized GenData structure which
 * 					on exit contains the training dataset
 * @param[in,out] 	test 		an initialized GenData structure which
 * 					on exit contains the test dataset
 * @param[in] 		cv_idx 		a vector of cv partitions created by
 * 					gensvm_make_cv_split()
 * @param[in] 		fold_idx 	index of the fold which becomes the
 * 					test dataset
 */
void gensvm_get_tt_split_dense(struct GenData *full,
		struct GenData *train, struct GenData *test,
		long *cv_idx, long fold_idx)
{
	long i, j, k, l, test_n, train_n;

	long n = full->n;
	long m = full->m;
	long K = full->K;

	double value;

	test_n = 0;
	for (i=0; i<n; i++)
		if (cv_idx[i] == fold_idx)
			test_n++;
	train_n = n - test_n;

	test->n = test_n;
	train->n = train_n;

	train->K = K;
	test->K = K;

	train->m = m;
	test->m = m;

	train->y = Calloc(long, train_n);
	test->y = Calloc(long, test_n);

	train->RAW = Calloc(double, train_n*(m+1));
	test->RAW = Calloc(double, test_n*(m+1));

	k = 0;
	l = 0;
	for (i=0; i<n; i++) {
		if (cv_idx[i] == fold_idx) {
			test->y[k] = full->y[i];
			for (j=0; j<m+1; j++) {
				value = matrix_get(full->RAW, n, m+1, i, j);
				matrix_set(test->RAW, test_n, m+1, k, j, 
						value);
			}
			k++;
		} else {
			train->y[l] = full->y[i];
			for (j=0; j<m+1; j++) {
				value = matrix_get(full->RAW, n, m+1, i, j);
				matrix_set(train->RAW, train_n, m+1, l, j, 
						value);
			}
			l++;
		}
	}

	train->Z = train->RAW;
	test->Z = test->RAW;
}


/**
 * @brief Create train and test dataset for a CV split with sparse data
 *
 * @details
 * Given a GenData structure for the full dataset, a previously created
 * cross validation split vector and a fold index, a training and test dataset
 * are created. It is assumed here that the data is stored as a sparse matrix,
 * and that the train and test data should also be stored as a sparse matrix.
 *
 * @sa
 * gensvm_get_tt_split_dense(), gensvm_get_tt_split()
 *
 * @param[in] 		full 		a GenData structure for the entire
 * 					dataset
 * @param[in,out] 	train 		an initialized GenData structure which
 * 					on exit contains the training dataset
 * @param[in,out] 	test 		an initialized GenData structure which
 * 					on exit contains the test dataset
 * @param[in] 		cv_idx 		a vector of cv partitions created by
 * 					gensvm_make_cv_split()
 * @param[in] 		fold_idx 	index of the fold which becomes the
 * 					test dataset
 */
void gensvm_get_tt_split_sparse(struct GenData *full,
		struct GenData *train, struct GenData *test,
		long *cv_idx, long fold_idx)
{
	long a, b, i, j, k, a_end, b_start, b_end, jx_elem, test_n, train_n, 
	     train_nnz, test_nnz, row_nnz, tr_nnz_idx, tr_row_idx, te_nnz_idx,
	     te_row_idx, tr_col_cnt, te_col_cnt;

	double value;
	enum {CSR, CSC} sp_type = full->spZ->type;

	// determine number of instances in test and train
	test_n = 0;
	for (i=0; i<full->n; i++)
		if (cv_idx[i] == fold_idx)
			test_n++;
	train_n = full->n - test_n;

	// set n, m, K variables
	train->n = train_n;
	train->m = full->m;
	train->K = full->K;
	test->n = test_n;
	test->m = full->m;
	test->K = full->K;

	// allocate outcome
	train->y = Calloc(long, train_n);
	test->y = Calloc(long, test_n);

	// compute train nnz and test nnz
	train_nnz = 0;
	test_nnz = 0;

	if (sp_type == CSR) {
		for (i=0; i<full->n; i++) {
			row_nnz = full->spZ->ix[i+1] - full->spZ->ix[i];
			if (cv_idx[i] == fold_idx) {
				test_nnz += row_nnz;
			} else {
				train_nnz += row_nnz;
			}
		}
	} else {
		for (k=0; k<full->spZ->nnz; k++) {
			i = full->spZ->jx[k];
			if (cv_idx[i] == fold_idx) {
				test_nnz++;
			} else {
				train_nnz++;
			}
		}
	}

	// allocate the train GenSparse
	train->spZ = gensvm_init_sparse();
	test->spZ = gensvm_init_sparse();

	// set the type the same as for the full dataset
	train->spZ->type = sp_type;
	test->spZ->type = sp_type;

	// set GenSparse variables for train
	train->spZ->nnz = train_nnz;
	train->spZ->n_row = train_n;
	train->spZ->n_col = full->spZ->n_col;
	train->spZ->values = Calloc(double, train_nnz);
	train->spZ->jx = Calloc(long, train_nnz);
	if (sp_type == CSR) {
		train->spZ->ix = Calloc(long, train_n+1);
	} else {
		train->spZ->ix = Calloc(long, train->spZ->n_col + 1);
	}

	// set GenSparse variables for test
	test->spZ->nnz = test_nnz;
	test->spZ->n_row = test_n;
	test->spZ->n_col = full->spZ->n_col;
	test->spZ->values = Calloc(double, test_nnz);
	test->spZ->jx = Calloc(long, test_nnz);
	if (sp_type == CSR) {
		test->spZ->ix = Calloc(long, test_n+1);
	} else {
		test->spZ->ix = Calloc(long, test->spZ->n_col + 1);
	}

	// copy over the labels
	j = 0;
	k = 0;
	for (i=0; i<full->n; i++) {
		if (cv_idx[i] == fold_idx)
			test->y[j++] = full->y[i];
		else
			train->y[k++] = full->y[i];
	}

	tr_nnz_idx = 0;
	tr_row_idx = 0;
	tr_col_cnt = 0;
	te_nnz_idx = 0;
	te_row_idx = 0;
	te_col_cnt = 0;

	test->spZ->ix[0] = 0;
	train->spZ->ix[0] = 0;

	a_end = (sp_type == CSR) ? full->spZ->n_row : 
		full->spZ->n_col;

	for (a=0; a<a_end; a++) {
		if (sp_type == CSC) {
			tr_col_cnt = 0;
			te_col_cnt = 0;
			tr_row_idx = 0;
			te_row_idx = 0;
		}

		b_start = full->spZ->ix[a];
		b_end = full->spZ->ix[a+1];

		for (b=b_start; b<b_end; b++) {
			value = full->spZ->values[b];
			i = (sp_type == CSR) ? a : full->spZ->jx[b];

			if (cv_idx[i] == fold_idx) {
				jx_elem = (sp_type == CSR) ? full->spZ->jx[b] : te_row_idx;
				test->spZ->values[te_nnz_idx] = value;
				test->spZ->jx[te_nnz_idx] = jx_elem;
				te_nnz_idx++;
				te_col_cnt++;
				if (sp_type == CSC)
					te_row_idx++;
			} else {
				jx_elem = (sp_type == CSR) ? full->spZ->jx[b] : tr_row_idx;
				train->spZ->values[tr_nnz_idx] = value;
				train->spZ->jx[tr_nnz_idx] = jx_elem;
				tr_nnz_idx++;
				tr_col_cnt++;
				if (sp_type == CSC)
					tr_row_idx++;
			}
		}

		if (sp_type == CSR) {
			if (cv_idx[a] == fold_idx) {
				test->spZ->ix[te_row_idx+1] = te_nnz_idx;
				te_row_idx++;
			} else {
				train->spZ->ix[tr_row_idx+1] = tr_nnz_idx;
				tr_row_idx++;
			}
		} else {
			test->spZ->ix[a+1] = test->spZ->ix[a] + te_col_cnt;
			train->spZ->ix[a+1] = train->spZ->ix[a] + tr_col_cnt;
		}
	}
}
