/**
 * @file globals.h
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Global definitions
 *
 * @details
 * This header file contains defines and includes which are used in many
 * parts of the program. Most notable are the Calloc, Malloc and Memset
 * defines, which are commonly used to allocate memory. These functions
 * are shorthands for their lowercase counterparts. 
 *
 * Furthermore, a maximum and minimum function are defined here. These 
 * functions have their own include guards, to ensure potential linked 
 * libraries don't conflict with these definitions.
 *
 */

#ifndef MSVMMAJ_GLOBALS_H
#define MSVMMAJ_GLOBALS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LENGTH 1024

#define Calloc(type, n) (type *)calloc((n), sizeof(type))
#define Malloc(type, n) (type *)malloc((n)*sizeof(type))
#define Memset(var, type, n) memset(var, 0, (n)*sizeof(type))

#ifndef MIN_MAX_DEFINE
#define MIN_MAX_DEFINE
#define maximum(a, b) (a) > (b) ? (a) : (b)
#define minimum(a, b) (a) < (b) ? (a) : (b)
#endif

#endif
