/**
 * @file gensvm_simplex.h
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Header file for gensvm_simplex.c
 *
 */

#ifndef GENSVM_SIMPLEX_H
#define GENSVM_SIMPLEX_H

// includes
#include "gensvm_base.h"

// forward declarations
void gensvm_simplex(struct GenModel *model);
void gensvm_simplex_diff(struct GenModel *model);

#endif
