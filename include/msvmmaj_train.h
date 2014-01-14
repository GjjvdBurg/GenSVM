/**
 * @file msvmmaj_train.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Header file for msvmmaj_train.c
 *
 * @details
 * Contains function declarations for functions used to train a single
 * MajModel.
 *
 */

#ifndef MSVMMAJ_TRAIN_H
#define MSVMMAJ_TRAIN_H

#include "globals.h"

//forward declarations
struct MajData;
struct MajModel;

// function declarations
void msvmmaj_optimize(struct MajModel *model, struct MajData *data);

double msvmmaj_get_loss(struct MajModel *model, struct MajData *data,
		double *ZV);

void msvmmaj_get_update(struct MajModel *model, struct MajData *data,
		double *B, double *ZAZ, double *ZAZV, double *ZAZVT);

#endif
