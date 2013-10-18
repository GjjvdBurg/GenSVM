#ifndef UTIL_H
#define UTIL_H

#include "globals.h"

// forward declarations
struct MajData;
struct MajModel;

// function declarations
void msvmmaj_read_data(struct MajData *dataset, char *data_file);

void msvmmaj_read_model(struct MajModel *model, char *model_filename);
void msvmmaj_write_model(struct MajModel *model, char *output_filename);

void msvmmaj_write_predictions(struct MajData *data, long *predy, 
		char *output_filename);

int msvmmaj_check_argv(int argc, char **argv, char *str);
int msvmmaj_check_argv_eq(int argc, char **argv, char *str);

void note(const char *fmt,...);

void msvmmaj_allocate_model(struct MajModel *model);
void msvmmaj_free_model(struct MajModel *model);
void msvmmaj_free_data(struct MajData *data);

#endif
