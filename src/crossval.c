#include "crossval.h"
#include "matrix.h"
#include "MSVMMaj.h"

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
						matrix_get(full_data->Z, m+1, i, j));
			k++;
		} else {
			train_data->y[l] = full_data->y[i];
			for (j=0; j<m+1; j++)
				matrix_set(train_data->Z, m+1, l, j,
						matrix_get(full_data->Z, m+1, i, j));
			l++;
		}
	}
}
