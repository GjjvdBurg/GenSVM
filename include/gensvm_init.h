/**
 * @file gensvm_init.h
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Header file for gensvm_init.c
 *
 * @details
 * Contains function declarations for the initialization functions for the
 * model weights and model V matrix.
 *
 * @copyright
 Copyright 2016, G.J.J. van den Burg.

 This file is part of GenSVM.

 GenSVM is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 GenSVM is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with GenSVM. If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef GENSVM_INIT_H
#define GENSVM_INIT_H

#include "gensvm_base.h"
#include "gensvm_rand.h"
#include "gensvm_print.h"

void gensvm_init_V(struct GenModel *from_model, struct GenModel *to_model,
		struct GenData *data);
void gensvm_initialize_weights(struct GenData *data, struct GenModel *model);

#endif
