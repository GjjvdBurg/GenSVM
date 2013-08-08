#include "util.h"
#include "matrix.h"
/*
	Read the data from the data_file. The data matrix X is augmented
	with a column of ones, to get the matrix Z.
*/
void read_data(struct Data *dataset, char *data_file)
{
	FILE *fid;
	long i, j;
	long n, m; // dimensions of data
	long nr = 0; // used to check consistency of data
	double value;
	long K = 0;
	long min_y = 1000;

	char buf[MAX_LINE_LENGTH];

	if ((fid = fopen(data_file, "r")) == NULL) {
		printf("\nERROR: datafile %s could not be opened.\n",
				data_file);
		exit(0);
	}

	// Read data dimensions
	nr += fscanf(fid, "%ld", &n);
	nr += fscanf(fid, "%ld", &m);

	// Allocate memory
	dataset->Z = Malloc(double, n*(m+1));

	// Read first line of data
	for (j=1; j<m+1; j++) {
		nr += fscanf(fid, "%lf", &value);
		matrix_set(dataset->Z, n, 0, j, value);
	}

	// Check if there is a label at the end of the line
	if (fgets(buf, MAX_LINE_LENGTH, fid) == NULL) {
		printf("ERROR: No label found on first line.\n");
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
			matrix_set(dataset->Z, m+1, i, j, value);
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
	} else if (min_y < 0 ) {
		printf("ERROR: wrong class labels in %s, minimum value is: %ld\n",
				data_file, min_y);
		exit(0);
	}

	if (nr < n * m) {
		printf("ERROR: not enough data found in %s\n", data_file);
		exit(0);
	}
	
	// Set the column of ones
	for (i=0; i<n; i++)
		matrix_set(dataset->Z, m+1, i, 0, 1.0);

	dataset->n = n;
	dataset->m = m;
	dataset->K = K;

	info("Succesfully read data file: %s\n",  data_file);
}

void next_line(FILE *fid, char *filename)
{
	char buffer[MAX_LINE_LENGTH];
	if (fgets(buffer, MAX_LINE_LENGTH, fid) == NULL) {
		fprintf(stderr, "Error reading file %s\n", filename);
		exit(1);
	}
}

double get_fmt_double(FILE *fid, char *filename, const char *fmt)
{
	char buffer[MAX_LINE_LENGTH];
	double value;

	if (fgets(buffer, MAX_LINE_LENGTH, fid) == NULL) {
		fprintf(stderr, "Error reading line from file %s\n", filename);
		exit(1);
	}
	sscanf(buffer, fmt, &value);

	return value;
}

long get_fmt_long(FILE *fid, char *filename, const char *fmt)
{
	char buffer[MAX_LINE_LENGTH];
	long value;

	if (fgets(buffer, MAX_LINE_LENGTH, fid) == NULL) {
		fprintf(stderr, "Error reading line from file %s\n", filename);
		exit(1);
	}
	sscanf(buffer, fmt, &value);

	return value;
}


void read_model(struct Model *model, char *model_filename)
{
	long i, j, nr = 0;
	FILE *fid;
	char buffer[MAX_LINE_LENGTH];
	char data_filename[MAX_LINE_LENGTH];
	double value = 0;

	fid = fopen(model_filename, "r");
	if (fid == NULL) {
		fprintf(stderr, "Error opening model file %s\n", model_filename);
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
	model->weight_idx = (int) get_fmt_long(fid, model_filename, "weight_idx = %li");

	// skip to data section
	for (i=0; i<2; i++)
		next_line(fid, model_filename);

	// read filename of data file
	if (fgets(buffer, MAX_LINE_LENGTH, fid) == NULL) {
		fprintf(stderr, "Error reading model file %s\n", model_filename);
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
				"Not enough elements of V found.\n", model_filename);
		exit(1);
	}

}

void write_model(struct Model *model, char *output_filename)
{
	FILE *fid;
	int i, j, diff, hours, minutes;
	char timestr[1000];
	time_t current_time, lt, gt;
	struct tm *lclt;

	// open output file
	fid = fopen(output_filename, "w");
	if (fid == NULL) {
		fprintf(stderr, "Error opening output file %s", output_filename);
		exit(1);
	}

	// get current time (in epoch)
	current_time = time(NULL);
	if (current_time == ((time_t)-1)) {
		fprintf(stderr, "Failed to compute the current time.\n");
		exit(1);
	}

	// convert time to local time and create a string
	lclt = localtime(&current_time);
	strftime(timestr, 1000, "%c", lclt);
	if (timestr == NULL) {
		fprintf(stderr, "Failed to convert time to string.\n");
		exit(1);
	}

	// calculate the difference from UTC including DST
	lt = mktime(localtime(&current_time));
	gt = mktime(gmtime(&current_time));
	diff = -difftime(gt, lt);
	hours = (diff/3600);
	minutes = (diff%3600)/60;
	if (lclt->tm_isdst == 1)
		hours++;

	// Write output to file
	fprintf(fid, "Output file for MSVMMaj (version %1.1f)\n", VERSION);
	fprintf(fid, "Generated on: %s (UTC %+03i:%02i)\n\n", timestr, hours, minutes);
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
			fprintf(fid, "%+15.16f ", matrix_get(model->V, model->K-1, i, j));
		}
		fprintf(fid, "\n");
	}

	fclose(fid);

}

void write_predictions(struct Data *data, long *predy, char *output_filename)
{
}

int check_argv(int argc, char **argv, char *str)
{
	int i;
	int arg_str = 0;
	for (i=1; i<argc; i++)
		if (strstr(argv[i], str) != NULL) {
			arg_str = i;
			break;
		}

	return arg_str;
}

int check_argv_eq(int argc, char **argv, char *str) 
{
	int i;
	int arg_str = 0;
	for (i=1; i<argc; i++)
		if (strcmp(argv[i], str) == 0) {
			arg_str = i;
			break;
		}

	return arg_str;
}

static void print_string_stdout(const char *s)
{
	fputs(s, stdout);
	fflush(stdout);
}

static void (*print_string) (const char *) = &print_string_stdout;

void set_print_string_function(void (*print_func)(const char *))
{
	if (print_func == NULL)
		print_string = &print_string_stdout;
	else
		print_string = print_func;
}

void info(const char *fmt,...)
{
	char buf[BUFSIZ];
	va_list ap;
	va_start(ap,fmt);
	vsprintf(buf,fmt,ap);
	va_end(ap);
	(*print_string)(buf);
}

double rnd()
{
	return (double) rand()/0x7FFFFFFF;
}

void allocate_model(struct Model *model)
{
	long n = model->n;
	long m = model->m;
	long K = model->K;

	model->W = Calloc(double, m*(K-1));
	if (model->W == NULL) {
		fprintf(stderr, "Failed to allocate memory for W.\n");
		exit(1);
	}

	model->t = Calloc(double, K-1);
	if (model->t == NULL) {
		fprintf(stderr, "Failed to allocate memory for t.\n");
		exit(1);
	}

	model->V = Calloc(double, (m+1)*(K-1));
	if (model->V == NULL) {
		fprintf(stderr, "Failed to allocate memory for V.\n");
		exit(1);
	}

	model->Vbar = Calloc(double, (m+1)*(K-1));
	if (model->Vbar == NULL) {
		fprintf(stderr, "Failed to allocate memory for Vbar.\n");
		exit(1);
	}

	model->U = Calloc(double, K*(K-1));
	if (model->U == NULL) {
		fprintf(stderr, "Failed to allocate memory for U.\n");
		exit(1);
	}

	model->UU = Calloc(double, n*K*(K-1));
	if (model->UU == NULL) {
		fprintf(stderr, "Failed to allocate memory for UU.\n");
		exit(1);
	}

	model->Q = Calloc(double, n*K);
	if (model->Q == NULL) {
		fprintf(stderr, "Failed to allocate memory for Q.\n");
		exit(1);
	}

	model->H = Calloc(double, n*K);
	if (model->H == NULL) {
		fprintf(stderr, "Failed to allocate memory for H.\n");
		exit(1);
	}

	model->R = Calloc(double, n*K);
	if (model->R == NULL) {
		fprintf(stderr, "Failed to allocate memory for R.\n");
		exit(1);
	}

	model->rho = Calloc(double, n);
	if (model->rho == NULL) {
		fprintf(stderr, "Failed to allocate memory for rho.\n");
		exit(1);
	}

}	

void free_model(struct Model *model)
{
	free(model->W);
	free(model->t);
	free(model->V);
	free(model->Vbar);
	free(model->U);
	free(model->UU);
	free(model->Q);
	free(model->H);
	free(model->rho);
	free(model->R);

	free(model);
}

void free_data(struct Data *data)
{
	free(data->Z);
	free(data->y);
	free(data);
}

/*
void print_matrix(double *M, long rows, long cols)
{
	long i, j;
	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			info("%8.8f ", matrix_get(M, cols, i, j));
		}
		info("\n");
	}
	info("\n");
}
*/
