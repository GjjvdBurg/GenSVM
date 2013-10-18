#ifndef GLOBALS_H
#define GLOBALS_H

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
