#include "util.h"

/*
	Read the data from the data_file. The data matrix X is augmented
	with a column of ones, to get the matrix Z.
*/
void read_data(struct Data *dataset, struct Model *model, char *data_file)
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

	model->n = n;
	model->m = m;
	model->K = K;

	info("Succesfully read data file: %s\n",  data_file);
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

/*
	Set a matrix element using ROW Major order. i denotes row,
	j denotes column.
*/
void matrix_set(double *M, long cols, long i, long j, double val)
{
	M[i*cols+j] = val;
}

/*
	Get a matrix element using ROW Major order. i denotes row,
	j denotes column.
*/
double matrix_get(double *M, long cols, long i, long j)
{
	return M[i*cols+j];
}

/*
	Add to an existing matrix element. Row-Major order is used.
*/
void matrix_add(double *M, long cols, long i, long j, double val)
{
	M[i*cols+j] += val;
}

/*
	Multiply existing matrix element. Row-Major order is used.
*/
void matrix_mult(double *M, long cols, long i, long j, double val)
{
	M[i*cols+j] *= val;
}


/*
	Set a matrix element of a 3D matrix in ROW major order. 
	N2 and N3 are the second and third dimension respectively
	and i, j, k are the indices of the first, second and third 
	dimensions respectively.
*/
void matrix3_set(double *M, long N2, long N3, long i, long j, long k, double val)
{
	M[k+N3*(j+N2*i)] = val;
}

/*
	Get a matrix element of a 3D matrix in ROW major order. 
	N2 and N3 are the second and third dimension respectively, and
	i, j and k are the indices of the first, second and third
	dimension of the requested element respectively.
*/
double matrix3_get(double *M, long N2, long N3, long i, long j, long k)
{
	return M[k+N3*(j+N2*i)];
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

void print_matrix(double *M, long rows, long cols)
{
	long i, j;
	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			printf("%8.8f ", matrix_get(M, cols, i, j));
		}
		printf("\n");
	}
	printf("\n");
}
