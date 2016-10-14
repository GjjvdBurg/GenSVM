/**
 * @file gensvm_optimize.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Header file for gensvm_optimize.c
 *
 * @details
 * Contains function declarations for functions used to train a single
 * GenModel.
 *
 */

#ifndef GENSVM_OPTIMIZE_H
#define GENSVM_OPTIMIZE_H

#include "gensvm_sv.h"
#include "gensvm_simplex.h"
#include "gensvm_update.h"

// function declarations
void gensvm_optimize(struct GenModel *model, struct GenData *data);
double gensvm_get_loss(struct GenModel *model, struct GenData *data, 
		struct GenWork *work);
void gensvm_calculate_errors(struct GenModel *model, struct GenData *data,
		double *ZV);
void gensvm_calculate_ZV_dense(struct GenModel *model, struct GenData *data,
		double *ZV);
void gensvm_calculate_ZV_sparse(struct GenModel *model, struct GenData *data,
		double *ZV);
void gensvm_calculate_huber(struct GenModel *model);
void gensvm_step_doubling(struct GenModel *model);

#endif
