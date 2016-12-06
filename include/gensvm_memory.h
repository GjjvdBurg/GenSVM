/**
 * @file gensvm_memory.h
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Header file for gensvm_memory.c
 *
 * @details
 * Contains macro definitions for easy memory access.

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

#ifndef GENSVM_MEMORY_H
#define GENSVM_MEMORY_H

#include <stddef.h>

/**
 * Wrapper macro for mycalloc(). This macro uses the __FILE__ and __LINE__ 
 * standard macros to fill in some of the arguments to mycalloc(). This macro 
 * should be used when writing code, not mycalloc().
 */
#define Calloc(type, size) \
	mycalloc(__FILE__, __LINE__, size, sizeof(type))

/**
 * Wrapper macro for mymalloc(). This macro uses the __FILE__ and __LINE__ 
 * standard macros to fill in some of the arguments to mymalloc(). This macro 
 * should be used when writing code, not mymalloc().
 */
#define Malloc(type, size) \
	mymalloc(__FILE__, __LINE__, (size)*sizeof(type))
/**
 * Wrapper macro for myrealloc(). This macro uses the __FILE__ and __LINE__ 
 * standard macros to fill in some of the arguments to myrealloc(). This macro 
 * should be used when writing code, not myrealloc().
 */
#define Realloc(var, type, size) \
	myrealloc(__FILE__, __LINE__, (size)*sizeof(type), var)
/**
 * Wrapper macro for memset(). Since memset is only used to zero a matrix, 
 * this macro is defined.
 */
#define Memset(var, type, size) \
	memset(var, 0, (size)*sizeof(type))

void *mycalloc(const char *file, int line, unsigned long size,
	size_t typesize);
void *mymalloc(const char *file, int line, unsigned long size);
void *myrealloc(const char *file, int line, unsigned long size, void *var);

#endif
