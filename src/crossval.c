/**
 * @file crossval.c
 * @author Gertjan van den Burg 
 * @date January 7, 2014
 * @brief Functions for cross validation
 *
 * @details
 * This file contains functions for performing cross validation. The funtion
 * msvmmaj_make_cv_split() creates a cross validation vector for non-stratified
 * cross validation. The function msvmmaj_get_tt_split() creates a train and 
 * test dataset from a given dataset and a pre-determined CV partition vector.
 * See individual function documentation for details.
 *
 */

#include "crossval.h"
#include "msvmmaj.h"
#include "msvmmaj_matrix.h"

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
void msvmmaj_make_cv_split(long N, long folds, long *cv_idx)
{
	long i, j, idx;

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
 * @brief Create train and test datasets for a CV split
 *
 * @details
 * Given a MajData structure for the full dataset, a previously created
 * cross validation split vector and a fold index, a training and test dataset
 * are created. 
 *
 * @param[in] 		full_data 	a MajData structure for the entire 
 * 					dataset
 * @param[in,out] 	train_data 	an initialized MajData structure which
 * 					on exit contains the training dataset
 * @param[in,out] 	test_data 	an initialized MajData structure which
 * 					on exit contains the test dataset
 * @param[in] 		cv_idx 		a vector of cv partitions created by
 * 					msvmmaj_make_cv_split()
 * @param[in] 		fold_idx 	index of the fold which becomes the 
 * 					test dataset
 */
void msvmmaj_get_tt_split(struct MajData *full_data, struct MajData *train_data,
		struct MajData *test_data, long *cv_idx, long fold_idx)
{
	long i, j, k, l, test_n, train_n;

	long n = full_data->n;
	long m = full_data->m;
	long K = full_data->K;

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

	train_data->Z = Calloc(double, train_n*(m+1));
	test_data->Z = Calloc(double, test_n*(m+1));

	k = 0;
	l = 0;
	for (i=0; i<n; i++) {
		if (cv_idx[i] == fold_idx) {
			test_data->y[k] = full_data->y[i];
			for (j=0; j<m+1; j++)
				matrix_set(test_data->Z, m+1, k, j, 
						matrix_get(full_data->Z, m+1, 
							i, j));
			k++;
		} else {
			train_data->y[l] = full_data->y[i];
			for (j=0; j<m+1; j++)
				matrix_set(train_data->Z, m+1, l, j,
						matrix_get(full_data->Z, m+1, 
							i, j));
			l++;
		}
	}
}
