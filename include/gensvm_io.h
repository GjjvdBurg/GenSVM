/**
 * @file gensvm_io.h
 * @author G.J.J. van den Burg
 * @date 2014-01-07
 * @brief Header file for gensvm_io.c
 *
 * @details
 * Function declarations for input/output functions.
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

#ifndef GENSVM_IO_H
#define GENSVM_IO_H

// includes
#include "gensvm_base.h"
#include "gensvm_print.h"
#include "gensvm_strutil.h"

// function declarations
void gensvm_read_data(struct GenData *dataset, char *data_file);
void gensvm_read_data_libsvm(struct GenData *dataset, char *data_file);

void gensvm_read_model(struct GenModel *model, char *model_filename);
void gensvm_write_model(struct GenModel *model, char *output_filename);

void gensvm_write_predictions(struct GenData *data, long *predy,
		char *output_filename);
void gensvm_time_string(char *buffer);

#endif
