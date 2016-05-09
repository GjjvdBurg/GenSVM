/**
 * @file gensvm_io.c
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Functions for input and output of data and model files
 *
 * @details
 * This file contains functions for reading and writing model files, and data
 * files. 
 *
 */

#include "gensvm.h"
#include "gensvm_io.h"
#include "gensvm_matrix.h"
#include "gensvm_strutil.h"
#include "gensvm_timer.h"

/**
 * @brief Read data from file
 * 
 * @details
 * Read the data from the data_file. The data matrix X is augmented
 * with a column of ones, to get the matrix Z. The data is expected
 * to follow a specific format, which is specified in the @ref spec_data_file.
 * The class labels are corrected internally to correspond to the interval 
 * [1 .. K], where K is the total number of classes.
 *
 * @todo
 * Make sure that this function allows datasets without class labels for
 * testing.
 *
 * @param[in,out] 	dataset 	initialized GenData struct
 * @param[in] 		data_file 	filename of the data file.
 */
void gensvm_read_data(struct GenData *dataset, char *data_file)
{
	FILE *fid;
	long i, j;
	long n, m; // dimensions of data
	long nr = 0; // used to check consistency of data
	double value;
	long K = 0;
	long min_y = 1000000;

	char buf[MAX_LINE_LENGTH];

	if ((fid = fopen(data_file, "r")) == NULL) {
		fprintf(stderr, "\nERROR: datafile %s could not be opened.\n",
				data_file);
		exit(0);
	}

	// Read data dimensions
	nr += fscanf(fid, "%ld", &n);
	nr += fscanf(fid, "%ld", &m);

	// Allocate memory
	dataset->RAW = Malloc(double, n*(m+1));

	// Read first line of data
	for (j=1; j<m+1; j++) {
		nr += fscanf(fid, "%lf", &value);
		matrix_set(dataset->RAW, n, 0, j, value);
	}

	// Check if there is a label at the end of the line
	if (fgets(buf, MAX_LINE_LENGTH, fid) == NULL) {
		fprintf(stderr, "ERROR: No label found on first line.\n");
		exit(1);
	}
	if (sscanf(buf, "%lf", &value) > 0) {
		dataset->y = Malloc(long, n);
		dataset->y[0] = value;
	} else if (dataset->y != NULL) {
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
			K = maximum(K, value);
			min_y = minimum(min_y, value);
		}
	}
	fclose(fid);

	// Correct labels: must be in [1, K]
	if (min_y == 0) {
		for (i=0; i<n; i++)
			dataset->y[i]++;
		K++;
	} else if (min_y < 0 ) {
		fprintf(stderr, "ERROR: wrong class labels in %s, minimum "
				"value is: %ld\n",
				data_file, min_y);
		exit(0);
	}

	if (nr < n * m) {
		fprintf(stderr, "ERROR: not enough data found in %s\n", 
				data_file);
		exit(0);
	}
	
	// Set the column of ones
	for (i=0; i<n; i++)
		matrix_set(dataset->RAW, m+1, i, 0, 1.0);

	dataset->n = n;
	dataset->m = m;
	dataset->r = m;
	dataset->K = K;
	dataset->Z = dataset->RAW;
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
	FILE *fid;
	char buffer[MAX_LINE_LENGTH];
	char data_filename[MAX_LINE_LENGTH];
	double value = 0;

	fid = fopen(model_filename, "r");
	if (fid == NULL) {
		fprintf(stderr, "Error opening model file %s\n", 
				model_filename);
		exit(1);
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
	if (fgets(buffer, MAX_LINE_LENGTH, fid) == NULL) {
		fprintf(stderr, "Error reading model file %s\n", 
				model_filename);
		exit(1);
	}
	sscanf(buffer, "filename = %s\n", data_filename);
	model->data_file = data_filename;

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
		fprintf(stderr, "Error reading model file %s. "
				"Not enough elements of V found.\n", 
				model_filename);
		exit(1);
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
	FILE *fid;
	long i, j;
	char timestr[MAX_LINE_LENGTH];

	// open output file
	fid = fopen(output_filename, "w");
	if (fid == NULL) {
		fprintf(stderr, "Error opening output file %s", 
				output_filename);
		exit(1);
	}
	get_time_string(timestr);

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
			fprintf(fid, "%+15.16f ", 
					matrix_get(model->V, 
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
	FILE *fid;

	fid = fopen(output_filename, "w");
	if (fid == NULL) {
		fprintf(stderr, "Error opening output file %s", 
				output_filename);
		exit(1);
	}
	
	fprintf(fid, "%li\n", data->n);
	fprintf(fid, "%li\n", data->m);

	for (i=0; i<data->n; i++) {
		for (j=0; j<data->m; j++) 
			fprintf(fid, "%f ", 
					matrix_get(data->Z, 
						data->m+1, i, j+1));
		fprintf(fid, "%li\n", predy[i]);
	}

	fclose(fid);
}
