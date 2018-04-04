/**
 * @file gensvm_print.h
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Header file for gensvm_print.c
 *
 * @details
 * Function declarations for printing to stdout and stderr.
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

#ifndef GENSVM_PRINT_H
#define GENSVM_PRINT_H

// includes
#include "gensvm_globals.h"

// function declarations
void note(const char *fmt,...);
void gensvm_error(const char *fmt,...);

// declare the function pointers used in note and err
void gensvm_print_output_fpt(const char *buf, ...);
void gensvm_print_error_fpt(const char *buf, ...);
void (*gensvm_print_out)(const char *, ...);
void (*gensvm_print_err)(const char *, ...);

#endif
