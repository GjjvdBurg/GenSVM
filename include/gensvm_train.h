/**
 * @file gensvm_train.h
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Header file for gensvm_train.c
 *
 */

#ifndef GENSVM_TRAIN_H
#define GENSVM_TRAIN_H

// includes
#include "gensvm_init.h"
#include "gensvm_kernel.h"
#include "gensvm_optimize.h"

// function declarations
void gensvm_train(struct GenModel *model, struct GenData *data,
		struct GenModel *seed_model);

#endif
