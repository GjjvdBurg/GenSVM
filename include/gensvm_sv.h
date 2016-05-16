/**
 * @file gensvm_sv.h
 * @author Gertjan van den Burg
 * @date May, 2014
 * @brief Header file for gensvm_sv.c
 *
 * @details
 * Contains function declarations for functions used to count support vectors.
 *
 */

#ifndef GENSVM_SV_H
#define GENSVM_SV_H

// includes
#include "gensvm_base.h"

// function declarations
long gensvm_num_sv(struct GenModel *model, struct GenData *data);

#endif
