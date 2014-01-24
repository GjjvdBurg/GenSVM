/**
 * @file msvmmaj_io.h
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Header files for msvmmaj_io.c
 *
 * @details
 * Function declarations for input/output functions.
 *
 */

#ifndef MSVMMAJ_IO_H
#define MSVMMAJ_IO_H

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

#endif
