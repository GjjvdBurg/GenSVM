/**
 * @file gensvm_cv_util.h
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Header file for gensvm_cv_util.c
 *
 * @details
 * Contains function declarations for functions needed for performing cross
 * validation on GenData structures.
 *
 */

#ifndef GENSVM_CV_UTIL_H
#define GENSVM_CV_UTIL_H

#include "gensvm_base.h"

void gensvm_make_cv_split(long N, long folds, long *cv_idx);
void gensvm_get_tt_split(struct GenData *full_data, struct GenData *train_data,
		struct GenData *test_data, long *cv_idx, long fold_idx);
void gensvm_get_tt_split_dense(struct GenData *full_data,
		struct GenData *train_data, struct GenData *test_data,
		long *cv_idx, long fold_idx);
void gensvm_get_tt_split_sparse(struct GenData *full_data,
		struct GenData *train_data, struct GenData *test_data,
		long *cv_idx, long fold_idx);

#endif
