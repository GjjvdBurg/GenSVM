/**
 * @file libGenSVM.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Header file for the core GenSVM library libGenSVM.c
 *
 * @details
 * The core computational routines for GenSVM are defined in libGenSVM.c.
 * This file contains function declarations for these functions.
 *
 */

/**
 * @todo
 * rename this file and libGenSVM.c to correspond with the lowercase convention.
 * Also change the name of the include guard.
 */
#ifndef LIBGENSVM_H
#define LIBGENSVM_H

#include "globals.h"

// forward declarations
struct GenData;
struct GenModel;

// function declarations
void gensvm_simplex_gen(long K, double *U);
void gensvm_category_matrix(struct GenModel *model, struct GenData *data);
void gensvm_simplex_diff(struct GenModel *model, struct GenData *dataset);

void gensvm_calculate_errors(struct GenModel *model, struct GenData *data,
	       	double *ZV);
void gensvm_calculate_huber(struct GenModel *model);

void gensvm_step_doubling(struct GenModel *model);

void gensvm_seed_model_V(struct GenModel *from_model,
	       	struct GenModel *to_model, struct GenData *data);
void gensvm_initialize_weights(struct GenData *data, struct GenModel *model);

#endif
