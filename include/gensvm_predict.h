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

#ifndef GENSVM_PREDICT_H
#define GENSVM_PREDICT_H

// includes
#include "gensvm_kernel.h"
#include "gensvm_simplex.h"
#include "gensvm_zv.h"

// function declarations
void gensvm_predict_labels(struct GenData *testdata,
	       	struct GenModel *model, long *predy);
void gensvm_predict_labels_dense(struct GenData *testdata,
		struct GenModel *model, long *predy);
void gensvm_predict_labels_sparse(struct GenData *testdata,
		struct GenModel *model, long *predy);
double gensvm_prediction_perf(struct GenData *data, long *perdy);

#endif
