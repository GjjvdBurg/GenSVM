/**
 * @file util.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Header file for util.c
 * 
 * @details
 * Function declarations for utility functions of the program.
 *
 */

#ifndef MSVMMAJ_UTIL_H
#define MSVMMAJ_UTIL_H

#include "globals.h"

// forward declarations
struct MajData;
struct MajModel;

// function declarations
int msvmmaj_check_argv(int argc, char **argv, char *str);
int msvmmaj_check_argv_eq(int argc, char **argv, char *str);

void note(const char *fmt,...);

#endif
