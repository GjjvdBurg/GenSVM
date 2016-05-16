/**
 * @file gensvm_init.h
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Header file for gensvm_init.c
 *
 * @details
 * Contains function declarations for the initialization functions for
 * GenModel and GenData structures.
 */

#ifndef GENSVM_INIT_H
#define GENSVM_INIT_H

// include
#include "globals.h"
#include "gensvm.h"

struct GenModel *gensvm_init_model();

struct GenData *gensvm_init_data();

void gensvm_allocate_model(struct GenModel *model);
void gensvm_reallocate_model(struct GenModel *model, long n, long m);
void gensvm_free_model(struct GenModel *model);
void gensvm_free_data(struct GenData *data);

#endif
