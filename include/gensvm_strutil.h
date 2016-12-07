/**
 * @file gensvm_strutil.h
 * @author G.J.J. van den Burg
 * @date 2013-08-01
 * @brief Header file for gensvm_strutil.c
 *
 * @details
 * Function declarations for useful string functions used in parsing
 * input files.
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

#ifndef GENSVM_STRUTIL_H
#define GENSVM_STRUTIL_H

#include "gensvm_globals.h"

bool str_startswith(const char *str, const char *pre);
bool str_endswith(const char *str, const char *suf);
bool str_contains_char(const char *str, const char c);
char **str_split(char *original, const char *delims, int *len_ret);

void next_line(FILE *fid, char *filename);
char *get_line(FILE *fid, char *filename, char *buffer);

double get_fmt_double(FILE *fid, char *filename, const char *fmt);
long get_fmt_long(FILE *fid, char *filename, const char *fmt);

long all_doubles_str(char *buffer, long offset, double *all_doubles);
long all_longs_str(char *buffer, long offset, long *all_longs);

#endif
