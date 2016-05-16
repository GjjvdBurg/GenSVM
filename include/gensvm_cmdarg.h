/**
 * @file gensvm_cmdarg.h
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Header file for gensvm_cmdarg.c
 *
 * @details
 * Function declarations for dealing with command line arguments.
 *
 */

#ifndef GENSVM_CMDARG_H
#define GENSVM_CMDARG_H

#include "globals.h"

// function declarations
int gensvm_check_argv(int argc, char **argv, char *str);
int gensvm_check_argv_eq(int argc, char **argv, char *str);

#endif
