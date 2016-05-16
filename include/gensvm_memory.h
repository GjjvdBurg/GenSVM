/**
 * @file gensvm_memory.h
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Global definitions
 *
 */

#ifndef GENSVM_MEMORY_H
#define GENSVM_MEMORY_H

#include <stddef.h>

#define Calloc(type, size) \
	mycalloc(__FILE__, __LINE__, size, sizeof(type))
#define Malloc(type, size) \
	mymalloc(__FILE__, __LINE__, (size)*sizeof(type))
#define Realloc(var, type, size) \
	myrealloc(__FILE__, __LINE__, (size)*sizeof(type), var)
#define Memset(var, type, size) \
	memset(var, 0, (size)*sizeof(type))

void *mycalloc(const char *file, int line, unsigned long size,
	size_t typesize);
void *mymalloc(const char *file, int line, unsigned long size);
void *myrealloc(const char *file, int line, unsigned long size, void *var);

#endif
