/**
 * @file gensvm_init.h
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Header file for gensvm_init.c
 *
 * @details
 * Contains function declarations for the initialization functions for the
 * model weights and model V matrix.
 */

#ifndef GENSVM_INIT_H
#define GENSVM_INIT_H

#include "gensvm_base.h"

void gensvm_init_V(struct GenModel *from_model, struct GenModel *to_model,
		struct GenData *data);
void gensvm_initialize_weights(struct GenData *data, struct GenModel *model);

#endif
