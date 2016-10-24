/**
 * @file gensvm_sv.h
 * @author G.J.J. van den Burg
 * @date 2014-05-01
 * @brief Header file for gensvm_sv.c
 *
 * @details
 * Contains function declarations for functions used to count support vectors.
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

#ifndef GENSVM_SV_H
#define GENSVM_SV_H

// includes
#include "gensvm_base.h"

// function declarations
long gensvm_num_sv(struct GenModel *model);

#endif
