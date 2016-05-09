/**
 * @file globals.h
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Global definitions
 *
 * @details
 * This header file contains defines and includes which are used in many
 * parts of the program. Most notably, it includes the gensvm_memory.h header
 * which defines functions for safe memory allocation.
 *
 * Furthermore, a maximum and minimum function are defined here. These
 * functions have their own include guards, to ensure potential linked
 * libraries don't conflict with these definitions.
 *
 */

#ifndef GENSVM_GLOBALS_H
#define GENSVM_GLOBALS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "gensvm_memory.h"

#define MAX_LINE_LENGTH 1024

#ifndef MIN_MAX_DEFINE
#define MIN_MAX_DEFINE
#define maximum(a, b) (a) > (b) ? (a) : (b)
#define minimum(a, b) (a) < (b) ? (a) : (b)
#endif

#endif
