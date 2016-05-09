/**
 * @file crossval.h
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Header file for crossval.c
 *
 * @details
 * Contains function declarations for functions needed for performing cross
 * validation on GenData structures.
 *
 */

#ifndef GENSVM_CROSSVAL_H
#define GENSVM_CROSSVAL_H

// forward delaration
struct GenData;

void gensvm_make_cv_split(long N, long folds, long *cv_idx);
void gensvm_get_tt_split(struct GenData *full_data, struct GenData *train_data,
		struct GenData *test_data, long *cv_idx, long fold_idx);

#endif
