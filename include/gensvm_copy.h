/**
 * @file gensvm_copy.h
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Header file for gensvm_copy.c
 *
 */

#ifndef GENSVM_COPY_H
#define GENSVM_COPY_H

// includes
#include "gensvm_base.h"

// function declarations
void gensvm_copy_model(struct GenModel *from, struct GenModel *to);

#endif
