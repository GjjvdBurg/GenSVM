#ifndef CROSSVAL_H
#define CROSSVAL_H

#include "globals.h"

// forward delaration
struct MajData;

void msvmmaj_make_cv_split(long N, long folds, long *cv_idx);
void msvmmaj_get_tt_split(struct MajData *full_data, struct MajData *train_data,
		struct MajData *test_data, long *cv_idx, long fold_idx);

#endif
