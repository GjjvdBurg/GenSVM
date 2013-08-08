#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "MSVMMaj.h"

#define Calloc(type, n) (type *)calloc((n), sizeof(type))
#define Malloc(type, n) (type *)malloc((n)*sizeof(type))
#define Memset(var, type, n) memset(var, 0, (n)*sizeof(type))
#define maximum(a, b) a > b ? a : b
#define minimum(a, b) a < b ? a : b

void read_data(struct Data *dataset, char *data_file);

void read_model(struct Model *model, char *model_filename);
void write_model(struct Model *model, char *output_filename);

void write_predictions(struct Data *data, long *predy, char *output_filename);

int check_argv(int argc, char **argv, char *str);
int check_argv_eq(int argc, char **argv, char *str);

void set_print_string_function(void (*print_func)(const char *));
void info(const char *fmt,...);

double rnd();

void allocate_model(struct Model *model);
void free_model(struct Model *model);
void free_data(struct Data *data);

