/**
 * @file strutil.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Header file for strutil.c
 *
 * @details
 * Function declarations for useful string functions used in parsing
 * input files.
 *
 */

#ifndef STRUTIL_H
#define STRUTIL_H

#include "globals.h"
#include "types.h"

bool str_startswith(const char *str, const char *pre);
bool str_endswith(const char *str, const char *suf);

void next_line(FILE *fid, char *filename);
void get_line(FILE *fid, char *filename, char *buffer);

double get_fmt_double(FILE *fid, char *filename, const char *fmt);
long get_fmt_long(FILE *fid, char *filename, const char *fmt);

long all_doubles_str(char *buffer, long offset, double *all_doubles);
long all_longs_str(char *buffer, long offset, long *all_longs);

#endif
