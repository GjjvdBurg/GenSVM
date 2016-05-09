/**
 * @file gensvm_io.h
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Header files for gensvm_io.c
 *
 * @details
 * Function declarations for input/output functions.
 *
 */

#ifndef GENSVM_IO_H
#define GENSVM_IO_H

// forward declarations
struct GenData;
struct GenModel;

// function declarations
void gensvm_read_data(struct GenData *dataset, char *data_file);

void gensvm_read_model(struct GenModel *model, char *model_filename);
void gensvm_write_model(struct GenModel *model, char *output_filename);

void gensvm_write_predictions(struct GenData *data, long *predy,
		char *output_filename);

#endif
