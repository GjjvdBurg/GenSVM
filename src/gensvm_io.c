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
		gensvm_error("[GenSVM Error]: Datafile %s could not be opened.\n",
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
		matrix_set(dataset->RAW, n, m+1, 0, j, value);
	}

	if (fgets(buf, GENSVM_MAX_LINE_LENGTH, fid) == NULL) {
		// LCOV_EXCL_START
		gensvm_error("[GenSVM Error]: No label found on first line.\n");
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
			matrix_set(dataset->RAW, n, m+1, i, j, value);
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
		gensvm_error("[GenSVM Error]: not enough data found in %s\n",
				data_file);
		exit(EXIT_FAILURE);
		// LCOV_EXCL_STOP
	}

	// Set the column of ones
	for (i=0; i<n; i++)
		matrix_set(dataset->RAW, n, m+1, i, 0, 1.0);

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
 * @brief Print an error to the screen and exit (copied from LibSVM)
 *
 * @param[in] 	line_num 	line number where the error occured
 *
 */
void exit_input_error(int line_num)
{
	gensvm_error("[GenSVM Error]: Wrong input format on line: %i\n", line_num);
	exit(EXIT_FAILURE);
}

/**
 * @brief Read data from a file in LibSVM/SVMlight format
 *
 * @details
 * This function reads data from a file where the data is stored in
 * LibSVM/SVMlight format. The file format is described in @ref 
 * spec_libsvm_data_file.  This is a sparse data format, which can be
 * beneficial for certain applications. The advantage of having this function
 * here is twofold: 1) existing datasets where data is stored in
 * LibSVM/SVMlight format can be easily used in GenSVM, and 2) sparse datasets
 * which are too large for memory when kept in dense format can be loaded
 * efficiently into GenSVM.
 *
 * @note
 * This code is based on the read_problem() function in the svm-train.c
 * file of LibSVM. It has however been expanded to be able to handle data 
 * files without labels.
 *
 * @note
 * This file tries to detect whether 1-based or 0-based indexing is used in 
 * the data file. By default 1-based indexing is used, but if an index is 
 * found with value 0, 0-based indexing is assumed.
 *
 * @sa
 * gensvm_read_problem()
 *
 * @param[in] 	 data 		GenData structure
 * @param[in] 	 data_file 	filename of the datafile
 *
 */
void gensvm_read_data_libsvm(struct GenData *data, char *data_file)
{
	bool do_sparse, zero_based = false;
	long i, j, n, m, K, nnz, cnt, tmp, index, row_cnt, num_labels, 
	     min_index = 1;
	int n_big, n_small, big_start;
	double value;
	FILE *fid = NULL;
	char *label = NULL,
	     *endptr = NULL,
	     **big_parts = NULL,
	     **small_parts = NULL;
	char buf[GENSVM_MAX_LINE_LENGTH];

	fid = fopen(data_file, "r");
	if (fid == NULL) {
		// LCOV_EXCL_START
		gensvm_error("[GenSVM Error]: Datafile %s could not be opened.\n",
				data_file);
		exit(EXIT_FAILURE);
		// LCOV_EXCL_STOP
	}

	// first count the number of elements
	n = 0;
	m = -1;

	num_labels = 0;
	nnz = 0;

	while (fgets(buf, GENSVM_MAX_LINE_LENGTH, fid) != NULL) {
		// split the string in labels and/or index:value pairs
		big_parts = str_split(buf, " \t", &n_big);

		// record if this line has a label (first part has no colon)
		num_labels += (!str_contains_char(big_parts[0], ':'));

		// check for each part if it is a index:value pair
		for (i=0; i<n_big; i++) {
			if (!str_contains_char(big_parts[i], ':'))
				continue;

			// split the index:value pair
			small_parts = str_split(big_parts[i], ":", &n_small);

			// convert the index to a number
			index = strtol(small_parts[0], &endptr, 10);

			// catch conversion errors
			if (endptr == small_parts[0] || errno != 0 ||
					*endptr != '\0')
				exit_input_error(n+1);

			// update the maximum index
			m = maximum(m, index);

			// update the minimum index
			min_index = minimum(min_index, index);

			// free the small parts
			for (j=0; j<n_small; j++) free(small_parts[j]);
			free(small_parts);

			// increment the nonzero counter
			nnz++;
		}

		// free the big parts
		for (i=0; i<n_big; i++) {
			free(big_parts[i]);
		}
		free(big_parts);

		// increment the number of observations
		n++;
	}

	// rewind the file pointer
	rewind(fid);

	// check if we have enough labels
	if (num_labels > 0 && num_labels != n) {
		gensvm_error("[GenSVM Error]: There are some lines with missing "
				"labels. Please fix this before "
				"continuing.\n");
		exit(EXIT_FAILURE);
	}

	// don't forget the column of ones
	nnz += n;

	// deal with 0-based or 1-based indexing in the LibSVM file
	if (min_index == 0) {
		m++;
		zero_based = true;
	}

	// check if sparsity is worth it
	do_sparse = gensvm_nnz_comparison(nnz, n, m+1);
	if (do_sparse) {
		data->spZ = gensvm_init_sparse();
		data->spZ->type = CSR;
		data->spZ->nnz = nnz;
		data->spZ->n_row = n;
		data->spZ->n_col = m+1;
		data->spZ->values = Calloc(double, nnz);

		if (data->spZ->type == CSR) {
			data->spZ->ix = Calloc(long, data->spZ->n_row+1);
		} else {
			data->spZ->ix = Calloc(long, data->spZ->n_col+1);
		}
		data->spZ->jx = Calloc(long, nnz);
		data->spZ->ix[0] = 0;
	} else {
		data->RAW = Calloc(double, n*(m+1));
		data->Z = data->RAW;
	}
	if (num_labels > 0)
		data->y = Calloc(long, n);

	K = 0;
	cnt = 0;
	for (i=0; i<n; i++) {
		if (fgets(buf, GENSVM_MAX_LINE_LENGTH, fid) == NULL) {
			// LCOV_EXCL_START
			gensvm_error("[GenSVM Error]: Error reading from data file %s\n",
					data_file);
			exit(EXIT_FAILURE);
			// LCOV_EXCL_STOP
		}

		// split the string in labels and/or index:value pairs
		big_parts = str_split(buf, " \t", &n_big);

		big_start = 0;
		// get the label from the first part if it exists
		if (!str_contains_char(big_parts[0], ':')) {
			label = strtok(big_parts[0], " \t\n");
			if (label == NULL) // empty line
				exit_input_error(i+1);

			// convert the label part to a number exit if there
			// are errors
			tmp = strtol(label, &endptr, 10);
			if (endptr == label || *endptr != '\0')
				exit_input_error(i+1);

			// assign label to y
			data->y[i] = tmp;

			// keep track of maximum K
			K = maximum(K, data->y[i]);

			// increment big part index
			big_start++;
		}

		row_cnt = 0;
		// set the first element in the row to 1 (this is the column 
		// of ones in the Z matrix for GenSVM)
		if (do_sparse) {
			data->spZ->values[cnt] = 1.0;
			data->spZ->jx[cnt] = 0;
			row_cnt++;
			cnt++;
		} else {
			matrix_set(data->RAW, n, m+1, i, 0, 1.0);
		}

		// read the rest of the line
		for (j=big_start; j<n_big; j++) {
			if (!str_contains_char(big_parts[j], ':'))
				continue;

			// split the index:value pair
			small_parts = str_split(big_parts[j], ":", &n_small);
			if (n_small != 2)
				exit_input_error(n+1);

			// convert the index to a long
			errno = 0;
			index = strtol(small_parts[0], &endptr, 10);

			// catch conversion errors
			if (endptr == small_parts[0] || errno != 0 ||
					*endptr != '\0')
				exit_input_error(n+1);

			// convert the value to a double
			errno = 0;
			value = strtod(small_parts[1], &endptr);
			if (endptr == small_parts[1] || errno != 0 ||
					(*endptr != '\0' && !isspace(*endptr)))
				exit_input_error(n+1);

			if (do_sparse) {
				data->spZ->values[cnt] = value;
				data->spZ->jx[cnt] = index + zero_based;
				row_cnt++;
				cnt++;
			} else {
				matrix_set(data->RAW, n, m+1, i, 
						index + zero_based, value);
			}

			// free the small parts
			free(small_parts[0]);
			free(small_parts[1]);
			free(small_parts);
		}

		if (do_sparse) {
			data->spZ->ix[i+1] = data->spZ->ix[i] + row_cnt;
		}

		// free the big parts
		for (j=0; j<n_big; j++) {
			free(big_parts[j]);
		}
		free(big_parts);
	}

	#if MAJOR_ORDER == 'c'
	if (do_sparse) {
		struct GenSparse *tmp_sparse = gensvm_sparse_csr_to_csc(
				data->spZ);
		gensvm_free_sparse(data->spZ);
		data->spZ = tmp_sparse;
	}
	#endif

	fclose(fid);

	data->n = n;
	data->m = m;
	data->r = m;
	data->K = K;

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
		gensvm_error("[GenSVM Error]: Couldn't open model file %s\n",
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
		gensvm_error("[GenSVM Error]: Error reading from model file %s\n",
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
			matrix_set(model->V, model->m+1, model->K-1, i, j, 
					value);
		}
	}
	if (nr != (model->m+1)*(model->K-1)) {
		// LCOV_EXCL_START
		gensvm_error("[GenSVM Error] Error reading from model file %s. "
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
		gensvm_error("[GenSVM Error]: Error opening output file %s\n",
				output_filename);
		exit(EXIT_FAILURE);
		// LCOV_EXCL_STOP
	}
	gensvm_time_string(timestr);

	// Write output to file
	fprintf(fid, "Output file for GenSVM (version %s)\n", VERSION_STRING);
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
						model->m+1, model->K-1, i, j));
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
		gensvm_error("[GenSVM Error]: Error opening output file %s\n",
				output_filename);
		exit(EXIT_FAILURE);
		// LCOV_EXCL_STOP
	}

	fprintf(fid, "%li\n", data->n);
	fprintf(fid, "%li\n", data->m);

	for (i=0; i<data->n; i++) {
		for (j=0; j<data->m; j++)
			fprintf(fid, "%.16f ", matrix_get(data->Z, data->n, 
						data->m+1, i, j+1));
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
		gensvm_error("[GenSVM Error]: Failed to compute the current time.\n");
		return;
		// LCOV_EXCL_STOP
	}

	// convert time to local time and create a string
	lclt = localtime(&current_time);
	size_t len = strftime(timestr, GENSVM_MAX_LINE_LENGTH, "%c", lclt);
	if (len == 0) {
		gensvm_error("[GenSVM Error]: Failed to convert time to string.\n");
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
