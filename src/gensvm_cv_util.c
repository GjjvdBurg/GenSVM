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
			idx = rand()%N;
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
 * This function tests if the data in the full_data structure is stored in a
 * dense matrix format or not, and calls gensvm_get_tt_split_dense() or
 * gensvm_get_tt_split_sparse() accordingly.
 *
 * @sa
 * gensvm_get_tt_split_dense(), gensvm_get_tt_split_sparse()
 *
 * @param[in] 		full_data 	a GenData structure for the entire
 * 					dataset
 * @param[in,out] 	train_data 	an initialized GenData structure which
 * 					on exit contains the training dataset
 * @param[in,out] 	test_data 	an initialized GenData structure which
 * 					on exit contains the test dataset
 * @param[in] 		cv_idx 		a vector of cv partitions created by
 * 					gensvm_make_cv_split()
 * @param[in] 		fold_idx 	index of the fold which becomes the
 * 					test dataset
 */
void gensvm_get_tt_split(struct GenData *full_data,
		struct GenData *train_data, struct GenData *test_data,
		long *cv_idx, long fold_idx)
{
	if (full_data->Z == NULL)
		gensvm_get_tt_split_sparse(full_data, train_data, test_data,
				cv_idx, fold_idx);
	else
		gensvm_get_tt_split_dense(full_data, train_data, test_data,
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
 * @param[in] 		full_data 	a GenData structure for the entire
 * 					dataset
 * @param[in,out] 	train_data 	an initialized GenData structure which
 * 					on exit contains the training dataset
 * @param[in,out] 	test_data 	an initialized GenData structure which
 * 					on exit contains the test dataset
 * @param[in] 		cv_idx 		a vector of cv partitions created by
 * 					gensvm_make_cv_split()
 * @param[in] 		fold_idx 	index of the fold which becomes the
 * 					test dataset
 */
void gensvm_get_tt_split_dense(struct GenData *full_data,
		struct GenData *train_data, struct GenData *test_data,
		long *cv_idx, long fold_idx)
{
	long i, j, k, l, test_n, train_n;

	long n = full_data->n;
	long m = full_data->m;
	long K = full_data->K;

	double value;

	test_n = 0;
	for (i=0; i<n; i++)
		if (cv_idx[i] == fold_idx)
			test_n++;
	train_n = n - test_n;

	test_data->n = test_n;
	train_data->n = train_n;

	train_data->K = K;
	test_data->K = K;

	train_data->m = m;
	test_data->m = m;

	train_data->y = Calloc(long, train_n);
	test_data->y = Calloc(long, test_n);

	train_data->RAW = Calloc(double, train_n*(m+1));
	test_data->RAW = Calloc(double, test_n*(m+1));

	k = 0;
	l = 0;
	for (i=0; i<n; i++) {
		if (cv_idx[i] == fold_idx) {
			test_data->y[k] = full_data->y[i];
			for (j=0; j<m+1; j++) {
				value = matrix_get(full_data->RAW, m+1, i, j);
				matrix_set(test_data->RAW, m+1, k, j, value);
			}
			k++;
		} else {
			train_data->y[l] = full_data->y[i];
			for (j=0; j<m+1; j++) {
				value = matrix_get(full_data->RAW, m+1, i, j);
				matrix_set(train_data->RAW, m+1, l, j, value);
			}
			l++;
		}
	}

	train_data->Z = train_data->RAW;
	test_data->Z = test_data->RAW;
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
 * @param[in] 		full_data 	a GenData structure for the entire
 * 					dataset
 * @param[in,out] 	train_data 	an initialized GenData structure which
 * 					on exit contains the training dataset
 * @param[in,out] 	test_data 	an initialized GenData structure which
 * 					on exit contains the test dataset
 * @param[in] 		cv_idx 		a vector of cv partitions created by
 * 					gensvm_make_cv_split()
 * @param[in] 		fold_idx 	index of the fold which becomes the
 * 					test dataset
 */
void gensvm_get_tt_split_sparse(struct GenData *full_data,
		struct GenData *train_data, struct GenData *test_data,
		long *cv_idx, long fold_idx)
{
	long i, j, test_n, train_n, train_nnz, test_nnz, row_nnz, jj,
	     jj_start, jj_end,
	     tr_nnz_idx = 0,
	     tr_row_idx = 0,
	     te_nnz_idx = 0,
	     te_row_idx = 0;

	double value;

	// determine number of instances in test and train
	test_n = 0;
	for (i=0; i<full_data->n; i++)
		if (cv_idx[i] == fold_idx)
			test_n++;
	train_n = full_data->n - test_n;

	// set n, m, K variables
	train_data->n = train_n;
	train_data->m = full_data->m;
	train_data->K = full_data->K;
	test_data->n = test_n;
	test_data->m = full_data->m;
	test_data->K = full_data->K;

	// allocate outcome
	train_data->y = Calloc(long, train_n);
	test_data->y = Calloc(long, test_n);

	// compute train nnz and test nnz
	train_nnz = 0;
	test_nnz = 0;
	for (i=0; i<full_data->n; i++) {
		row_nnz = full_data->spZ->ia[i+1] - full_data->spZ->ia[i];
		if (cv_idx[i] == fold_idx) {
			test_nnz += row_nnz;
		} else {
			train_nnz += row_nnz;
		}
	}

	// allocate the train GenSparse
	train_data->spZ = gensvm_init_sparse();
	test_data->spZ = gensvm_init_sparse();

	// set GenSparse variables for train
	train_data->spZ->nnz = train_nnz;
	train_data->spZ->n_row = train_n;
	train_data->spZ->n_col = full_data->m+1;
	train_data->spZ->values = Calloc(double, train_nnz);
	train_data->spZ->ia = Calloc(long, train_n+1);
	train_data->spZ->ja = Calloc(long, train_nnz);

	// set GenSparse variables for test
	test_data->spZ->nnz = test_nnz;
	test_data->spZ->n_row = test_n;
	test_data->spZ->n_col = full_data->m+1;
	test_data->spZ->values = Calloc(double, test_nnz);
	test_data->spZ->ia = Calloc(long, test_n+1);
	test_data->spZ->ja = Calloc(long, test_nnz);

	tr_nnz_idx = 0;
	tr_row_idx = 0;
	te_nnz_idx = 0;
	te_row_idx = 0;

	test_data->spZ->ia[0] = 0;
	train_data->spZ->ia[0] = 0;
	for (i=0; i<full_data->n; i++) {
		jj_start = full_data->spZ->ia[i];
		jj_end = full_data->spZ->ia[i+1];

		for (jj=jj_start; jj<jj_end; jj++) {
			j = full_data->spZ->ja[jj];
			value = full_data->spZ->values[jj];

			if (cv_idx[i] == fold_idx) {
				test_data->spZ->values[te_nnz_idx] = value;
				test_data->spZ->ja[te_nnz_idx] = j;
				te_nnz_idx++;
			} else {
				train_data->spZ->values[tr_nnz_idx] = value;
				train_data->spZ->ja[tr_nnz_idx] = j;
				tr_nnz_idx++;
			}
		}

		if (cv_idx[i] == fold_idx) {
			test_data->y[te_row_idx] = full_data->y[i];
			test_data->spZ->ia[te_row_idx+1] = te_nnz_idx;
			te_row_idx++;
		} else {
			train_data->y[tr_row_idx] = full_data->y[i];
			train_data->spZ->ia[tr_row_idx+1] = tr_nnz_idx;
			tr_row_idx++;
		}
	}
}
