/**
 * @file gensvm_pred.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Header file for gensvm_pred.c
 *
 * @details
 * Contains function declarations for prediction functions.
 *
 */

#ifndef GENSVM_PRED_H
#define GENSVM_PRED_H

#include "globals.h"

// forward declarations
struct GenData;
struct GenModel;

// function declarations
void gensvm_predict_labels(struct GenData *testdata,
	       	struct GenModel *model, long *predy);
double gensvm_prediction_perf(struct GenData *data, long *perdy);

#endif
