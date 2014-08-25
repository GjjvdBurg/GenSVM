/**
 * @file gensvm_train.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Header file for gensvm_train.c
 *
 * @details
 * Contains function declarations for functions used to train a single
 * GenModel.
 *
 */

#ifndef GENSVM_TRAIN_H
#define GENSVM_TRAIN_H

#include "globals.h"

//forward declarations
struct GenData;
struct GenModel;

// function declarations
void gensvm_optimize(struct GenModel *model, struct GenData *data);

double gensvm_get_loss(struct GenModel *model, struct GenData *data,
		double *ZV);

void gensvm_get_update(struct GenModel *model, struct GenData *data,
		double *B, double *ZAZ, double *ZAZV, double *ZAZVT);

#endif
