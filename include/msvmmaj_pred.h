/**
 * @file msvmmaj_pred.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Header file for msvmmaj_pred.c
 *
 * @details
 * Contains function declarations for prediction functions.
 *
 */

#ifndef MSVMMAJ_PRED_H
#define MSVMMAJ_PRED_H

#include "globals.h"

// forward declarations
struct MajData;
struct MajModel;

// function declarations
void msvmmaj_predict_labels(struct MajData *data, struct MajModel *model,
		long *predy);
double msvmmaj_prediction_perf(struct MajData *data, long *perdy);

#endif
