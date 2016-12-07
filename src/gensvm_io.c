/**
 * @file gensvm_io.c
 * @author G.J.J. van den Burg
 * @date 2014-01-07
 * @brief Functions for input and output of data and model files
 *
 * @details
 * This file contains functions for reading and writing model files, and data
 * files. It also contains a function for generating a string of the current 
 * time, used in writing output files.
 *
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

#include "gensvm_io.h"

/**
 * @brief Read data from file
 *
 * @details
 * Read the data from the data_file. The data matrix X is augmented
 * with a column of ones, to get the matrix Z. The data is expected
 * to follow a specific format, which is specified in the @ref spec_data_file.
 * The class labels are assumed to be in the interval [1 .. K], which can be 
 * checked using the function gensvm_check_outcome_contiguous().
 *
 * @param[in,out] 	dataset 	initialized GenData struct
 * @param[in] 		data_file 	filename of the data file.
 */
void gensvm_read_data(struct GenData *dataset, char *data_file)
{
	FILE *fid = NULL;
	long i, j, n, m,
	     nr = 0,
	     K = 0;
	double value;
	char buf[GENSVM_MAX_LINE_LENGTH];

	if ((fid = fopen(data_file, "r")) == NULL) {
		// LCOV_EXCL_START
		err("[GenSVM Error]: Datafile %s could not be opened.\n",
				data_file);
		exit(EXIT_FAILURE);
		// LCOV_EXCL_STOP
	}

	// Read data dimensions
	nr += fscanf(fid, "%ld", &n);
	nr += fscanf(fid, "%ld", &m);

	// Allocate memory
	dataset->RAW = Malloc(double, n*(m+1));

	// Read first line of data
	for (j=1; j<m+1; j++) {
		nr += fscanf(fid, "%lf", &value);
		matrix_set(dataset->RAW, m+1, 0, j, value);
	}

	if (fgets(buf, GENSVM_MAX_LINE_LENGTH, fid) == NULL) {
		// LCOV_EXCL_START
		err("[GenSVM Error]: No label found on first line.\n");
		exit(EXIT_FAILURE);
		// LCOV_EXCL_STOP
	}

	// Check if there is a label at the end of the line
	if (sscanf(buf, "%lf", &value) > 0) {
		dataset->y = Malloc(long, n);
		dataset->y[0] = value;
		K = 1;
	} else {
		free(dataset->y);
		dataset->y = NULL;
	}

	// Read the rest of the file
	for (i=1; i<n; i++) {
		for (j=1; j<m+1; j++) {
			nr += fscanf(fid, "%lf", &value);
			matrix_set(dataset->RAW, m+1, i, j, value);
		}
		if (dataset->y != NULL) {
			nr += fscanf(fid, "%lf", &value);
			dataset->y[i] = (long) value;
			K = maximum(K, dataset->y[i]);
		}
	}
	fclose(fid);

	if (nr < n * m) {
		// LCOV_EXCL_START
		err("[GenSVM Error]: not enough data found in %s\n",
				data_file);
		exit(EXIT_FAILURE);
		// LCOV_EXCL_STOP
	}

	// Set the column of ones
	for (i=0; i<n; i++)
		matrix_set(dataset->RAW, m+1, i, 0, 1.0);

	dataset->n = n;
	dataset->m = m;
	dataset->r = m;
	dataset->K = K;
	dataset->Z = dataset->RAW;

	if (gensvm_could_sparse(dataset->Z, n, m+1)) {
		note("Converting to sparse ... ");
		dataset->spZ = gensvm_dense_to_sparse(dataset->Z, n, m+1);
		note("done.\n");
		free(dataset->RAW);
		dataset->RAW = NULL;
		dataset->Z = NULL;
	}
}


/**
 * @brief Read model from file
 *
 * @details
 * Read a GenModel from a model file. The GenModel struct must have been
 * initalized elswhere. The model file is expected to follow the @ref
 * spec_model_file. The easiest way to generate a model file is through
 * gensvm_write_model(), which can for instance be used in trainGenSVM.c.
 *
 * @param[in,out] 	model 		initialized GenModel
 * @param[in] 		model_filename 	filename of the model file
 *
 */
void gensvm_read_model(struct GenModel *model, char *model_filename)
{
	long i, j, nr = 0;
	FILE *fid = NULL;
	char buffer[GENSVM_MAX_LINE_LENGTH];
	char data_filename[GENSVM_MAX_LINE_LENGTH];
	double value = 0;

	fid = fopen(model_filename, "r");
	if (fid == NULL) {
		// LCOV_EXCL_START
		err("[GenSVM Error]: Couldn't open model file %s\n",
				model_filename);
		exit(EXIT_FAILURE);
		// LCOV_EXCL_STOP
	}
	// skip the first four lines
	for (i=0; i<4; i++)
		next_line(fid, model_filename);

	// read all model variables
	model->p = get_fmt_double(fid, model_filename, "p = %lf");
	model->lambda = get_fmt_double(fid, model_filename, "lambda = %lf");
	model->kappa = get_fmt_double(fid, model_filename, "kappa = %lf");
	model->epsilon = get_fmt_double(fid, model_filename, "epsilon = %lf");
	model->weight_idx = (int) get_fmt_long(fid, model_filename,
			"weight_idx = %li");

	// skip to data section
	for (i=0; i<2; i++)
		next_line(fid, model_filename);

	// read filename of data file
	if (fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid) == NULL) {
		// LCOV_EXCL_START
		err("[GenSVM Error]: Error reading from model file %s\n",
				model_filename);
		exit(EXIT_FAILURE);
		// LCOV_EXCL_STOP
	}
	sscanf(buffer, "filename = %s\n", data_filename);
	model->data_file = Calloc(char, GENSVM_MAX_LINE_LENGTH);
	strcpy(model->data_file, data_filename);

	// read all data variables
	model->n = get_fmt_long(fid, model_filename, "n = %li\n");
	model->m = get_fmt_long(fid, model_filename, "m = %li\n");
	model->K = get_fmt_long(fid, model_filename, "K = %li\n");

	// skip to output
	for (i=0; i<2; i++)
		next_line(fid, model_filename);

	// read the matrix V and check for consistency
	model->V = Malloc(double, (model->m+1)*(model->K-1));
	for (i=0; i<model->m+1; i++) {
		for (j=0; j<model->K-1; j++) {
			nr += fscanf(fid, "%lf ", &value);
			matrix_set(model->V, model->K-1, i, j, value);
		}
	}
	if (nr != (model->m+1)*(model->K-1)) {
		// LCOV_EXCL_START
		err("[GenSVM Error] Error reading from model file %s. "
				"Not enough elements of V found.\n",
				model_filename);
		exit(EXIT_FAILURE);
		// LCOV_EXCL_STOP
	}
}

/**
 * @brief Write model to file
 *
 * @details
 * Write a GenModel to a file. The current time is specified in the file in
 * UTC + offset. The model file further corresponds to the @ref
 * spec_model_file.
 *
 * @param[in] 	model 		GenModel which contains an estimate for
 * 				GenModel::V
 * @param[in] 	output_filename the output file to write the model to
 *
 */
void gensvm_write_model(struct GenModel *model, char *output_filename)
{
	FILE *fid = NULL;
	long i, j;
	char timestr[GENSVM_MAX_LINE_LENGTH];

	// open output file
	fid = fopen(output_filename, "w");
	if (fid == NULL) {
		// LCOV_EXCL_START
		err("[GenSVM Error]: Error opening output file %s\n",
				output_filename);
		exit(EXIT_FAILURE);
		// LCOV_EXCL_STOP
	}
	gensvm_time_string(timestr);

	// Write output to file
	fprintf(fid, "Output file for GenSVM (version %1.1f)\n", VERSION);
	fprintf(fid, "Generated on: %s\n\n", timestr);
	fprintf(fid, "Model:\n");
	fprintf(fid, "p = %15.16f\n", model->p);
	fprintf(fid, "lambda = %15.16f\n", model->lambda);
	fprintf(fid, "kappa = %15.16f\n", model->kappa);
	fprintf(fid, "epsilon = %g\n", model->epsilon);
	fprintf(fid, "weight_idx = %i\n", model->weight_idx);
	fprintf(fid, "\n");
	fprintf(fid, "Data:\n");
	fprintf(fid, "filename = %s\n", model->data_file);
	fprintf(fid, "n = %li\n", model->n);
	fprintf(fid, "m = %li\n", model->m);
	fprintf(fid, "K = %li\n", model->K);
	fprintf(fid, "\n");
	fprintf(fid, "Output:\n");
	for (i=0; i<model->m+1; i++) {
		for (j=0; j<model->K-1; j++) {
			if (j > 0)
				fprintf(fid, " ");
			fprintf(fid, "%+15.16f", matrix_get(model->V,
						model->K-1, i, j));
		}
		fprintf(fid, "\n");
	}

	fclose(fid);
}

/**
 * @brief Write predictions to file
 *
 * @details
 * Write the given predictions to an output file, such that the resulting file
 * corresponds to the @ref spec_data_file.
 *
 * @param[in] 	data 		GenData with the original instances
 * @param[in] 	predy 		predictions of the class labels of the
 * 				instances in the given GenData. Note that the
 * 				order of the instances is assumed to be the
 * 				same.
 * @param[in] 	output_filename the file to which the predictions are written
 *
 */
void gensvm_write_predictions(struct GenData *data, long *predy,
		char *output_filename)
{
	long i, j;
	FILE *fid = NULL;

	fid = fopen(output_filename, "w");
	if (fid == NULL) {
		// LCOV_EXCL_START
		err("[GenSVM Error]: Error opening output file %s\n",
				output_filename);
		exit(EXIT_FAILURE);
		// LCOV_EXCL_STOP
	}

	fprintf(fid, "%li\n", data->n);
	fprintf(fid, "%li\n", data->m);

	for (i=0; i<data->n; i++) {
		for (j=0; j<data->m; j++)
			fprintf(fid, "%.16f ", matrix_get(data->Z, data->m+1, i,
						j+1));
		fprintf(fid, "%li\n", predy[i]);
	}

	fclose(fid);
}

/**
 * @brief Get time string with UTC offset
 *
 * @details
 * Create a string for the current system time. Include an offset of UTC for
 * consistency. The format of the generated string is "DDD MMM D HH:MM:SS
 * YYYY (UTC +HH:MM)", e.g. "Fri Aug 9, 12:34:56 2013 (UTC +02:00)".
 *
 * @param[in,out] 	buffer 	allocated string buffer, on exit contains
 * 				formatted string
 *
 */
void gensvm_time_string(char *buffer)
{
	int diff, hours, minutes;
	char timestr[GENSVM_MAX_LINE_LENGTH];
	time_t current_time, lt, gt;
	struct tm *lclt = NULL;

	// get current time (in epoch)
	current_time = time(NULL);
	if (current_time == ((time_t)-1)) {
		// LCOV_EXCL_START
		err("[GenSVM Error]: Failed to compute the current time.\n");
		return;
		// LCOV_EXCL_STOP
	}

	// convert time to local time and create a string
	lclt = localtime(&current_time);
	strftime(timestr, GENSVM_MAX_LINE_LENGTH, "%c", lclt);
	if (timestr == NULL) {
		err("[GenSVM Error]: Failed to convert time to string.\n");
		return;
	}

	// calculate the UTC offset including DST
	lt = mktime(localtime(&current_time));
	gt = mktime(gmtime(&current_time));
	diff = -difftime(gt, lt);
	hours = (diff/3600);
	minutes = (diff%3600)/60;
	if (lclt->tm_isdst == 1)
		hours++;

	sprintf(buffer, "%s (UTC %+03i:%02i)", timestr, hours, minutes);
}
