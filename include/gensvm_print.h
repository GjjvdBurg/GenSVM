/**
 * @file gensvm_print.h
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Header file for gensvm_print.c
 *
 * @details
 * Function declarations for printing to stdout and stderr.
 *
 */

#ifndef GENSVM_PRINT_H
#define GENSVM_PRINT_H

// includes
#include "globals.h"

// function declarations
void note(const char *fmt,...);
void err(const char *fmt,...);

#endif
