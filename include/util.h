#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "MSVMMaj.h"

#define Calloc(type, n) (type *)calloc((n), sizeof(type))
#define Malloc(type, n) (type *)malloc((n)*sizeof(type))
#define Memset(var, type, n) memset(var, 0, (n)*sizeof(type))
#define maximum(a, b) a > b ? a : b
#define minimum(a, b) a < b ? a : b

void read_data(struct Data *dataset, struct Model *model, char *data_file);

int check_argv(int argc, char **argv, char *str);
int check_argv_eq(int argc, char **argv, char *str);

void set_print_string_function(void (*print_func)(const char *));
void info(const char *fmt,...);

double rnd();

void matrix_set(double *M, long cols, long i, long j, double val);
void matrix_add(double *M, long cols, long i, long j, double val);
void matrix_mult(double *M, long cols, long i, long j, double val);
double matrix_get(double *M, long cols, long i, long j);

void matrix3_set(double *M, long N2, long N3, long i, long j, long k, double val);
double matrix3_get(double *M, long N2, long N3, long i, long j, long k);

void allocate_model(struct Model *model);
void free_model(struct Model *model);
void free_data(struct Data *data);

void print_matrix(double *M, long rows, long cols);
