/**
 * @file GenSVM_grid.c
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Command line interface for the grid search program
 *
 * @details
 * This is a command line interface to the parameter grid search functionality
 * of the algorithm. The grid search is specified in a separate file, thereby
 * reducing the number of command line arguments. See
 * read_grid_from_file() for documentation on the grid file.
 *
 * The program runs a grid search as specified in the grid file. If
 * desired the grid search can incorporate consistency checks to find the
 * configuration among the best configurations which scores consistently high.
 * All output is written to stdout, unless the quiet mode is specified.
 *
 * For further usage information, see the program help function.
 *
 */

#include "gensvm_cmdarg.h"
#include "gensvm_io.h"
#include "gensvm_gridsearch.h"

#define MINARGS 2

extern FILE *GENSVM_OUTPUT_FILE;
extern FILE *GENSVM_ERROR_FILE;

// function declarations
void exit_with_help();
void parse_command_line(int argc, char **argv, char *input_filename);
void read_grid_from_file(char *input_filename, struct GenGrid *grid);

/**
 * @brief Help function
 */
void exit_with_help()
{
	printf("This is GenSVM, version %1.1f\n\n", VERSION);
	printf("Usage: trainGenSVMdataset [options] grid_file\n");
	printf("Options:\n");
	printf("-h | -help : print this help.\n");
	printf("-q : quiet mode (no output, not even errors!)\n");

	exit(EXIT_FAILURE);
}

/**
 * @brief Main interface function for trainGenSVMdataset
 *
 * @details
 * Main interface for the command line program. A given grid file which
 * specifies a grid search over a single dataset is read. From this, a Queue
 * is created containing all Task instances that need to be performed in the
 * search. Depending on the type of dataset, either cross validation or
 * train/test split training is performed for all tasks. If specified,
 * consistency repeats are done at the end of the grid search. Note that
 * currently no output is produced other than what is written to stdout.
 *
 * @param[in] 	argc 	number of command line arguments
 * @param[in] 	argv 	array of command line arguments
 *
 */
int main(int argc, char **argv)
{
	char input_filename[MAX_LINE_LENGTH];

	struct GenGrid *grid = gensvm_init_grid();
	struct GenData *train_data = gensvm_init_data();
	struct GenData *test_data = gensvm_init_data();
	struct GenQueue *q = gensvm_init_queue();

	if (argc < MINARGS || gensvm_check_argv(argc, argv, "-help")
			|| gensvm_check_argv_eq(argc, argv, "-h") )
		exit_with_help();
	parse_command_line(argc, argv, input_filename);

	note("Reading grid file\n");
	read_grid_from_file(input_filename, grid);

	note("Reading data from %s\n", grid->train_data_file);
	gensvm_read_data(train_data, grid->train_data_file);
	if (grid->traintype == TT) {
		note("Reading data from %s\n", grid->test_data_file);
		gensvm_read_data(test_data, grid->test_data_file);
	}

	note("Creating queue\n");
	gensvm_fill_queue(grid, q, train_data, test_data);

	srand(time(NULL));

	note("Starting training\n");
	start_training(q);
	note("Training finished\n");

	if (grid->repeats > 0) {
		consistency_repeats(q, grid->repeats, grid->traintype);
	}

	gensvm_free_queue(q);
	gensvm_free_grid(grid);
	gensvm_free_data(train_data);
	gensvm_free_data(test_data);

	note("Done.\n");
	return 0;
}

/**
 * @brief Parse command line arguments
 *
 * @details
 * Few arguments can be supplied to the command line. Only quiet mode can be
 * specified, or help can be requested. The filename of the grid file is
 * read from the arguments. Parsing of the grid file is done separately in
 * read_grid_from_file().
 *
 * @param[in] 	argc 		number of command line arguments
 * @param[in] 	argv 		array of command line arguments
 * @param[in] 	input_filename 	pre-allocated buffer for the grid
 * 				filename.
 *
 */
void parse_command_line(int argc, char **argv, char *input_filename)
{
	int i;

	GENSVM_OUTPUT_FILE = stdout;
	GENSVM_ERROR_FILE = stderr;

	for (i=1; i<argc; i++) {
		if (argv[i][0] != '-') break;
		if (++i>=argc)
			exit_with_help();
		switch (argv[i-1][1]) {
			case 'q':
				GENSVM_OUTPUT_FILE = NULL;
				GENSVM_ERROR_FILE = NULL;
				i--;
				break;
			default:
				fprintf(stderr, "Unknown option: -%c\n",
						argv[i-1][1]);
				exit_with_help();
		}
	}

	if (i >= argc)
		exit_with_help();

	strcpy(input_filename, argv[i]);
}

KernelType parse_kernel_str(char *kernel_line)
{
	if (str_endswith(kernel_line, "LINEAR\n")) {
		return K_LINEAR;
	} else if (str_endswith(kernel_line, "POLY\n")) {
		return K_POLY;
	} else if (str_endswith(kernel_line, "RBF\n")) {
		return K_RBF;
	} else if (str_endswith(kernel_line, "SIGMOID\n")) {
		return K_SIGMOID;
	} else {
		fprintf(stderr, "Unknown kernel specified on line: %s\n",
				kernel_line);
		exit(1);
	}
}

/**
 * @brief Read the GenGrid struct from file
 *
 * @details
 * Read the GenGrid struct from a file. The grid file follows a specific
 * format specified in @ref spec_grid_file.
 *
 * Commonly used string functions in this function are all_doubles_str() and
 * all_longs_str().
 *
 * @param[in] 	input_filename 	filename of the grid file
 * @param[in] 	grid 	GenGrid structure to place the parsed
 * 				parameter grid.
 *
 */
void read_grid_from_file(char *input_filename, struct GenGrid *grid)
{
	long i, nr = 0;
	FILE *fid;
	char buffer[MAX_LINE_LENGTH];
	char train_filename[MAX_LINE_LENGTH];
	char test_filename[MAX_LINE_LENGTH];
	double *params = Calloc(double, MAX_LINE_LENGTH);
	long *lparams = Calloc(long, MAX_LINE_LENGTH);

	fid = fopen(input_filename, "r");
	if (fid == NULL) {
		fprintf(stderr, "Error opening grid file %s\n",
				input_filename);
		exit(1);
	}
	grid->traintype = CV;
	while ( fgets(buffer, MAX_LINE_LENGTH, fid) != NULL ) {
		Memset(params, double,  MAX_LINE_LENGTH);
		Memset(lparams, long, MAX_LINE_LENGTH);
		if (str_startswith(buffer, "train:")) {
			sscanf(buffer, "train: %s\n", train_filename);
			grid->train_data_file = Calloc(char,
					MAX_LINE_LENGTH);
			strcpy(grid->train_data_file, train_filename);
		} else if (str_startswith(buffer, "test:")) {
			sscanf(buffer, "test: %s\n", test_filename);
			grid->test_data_file = Calloc(char,
				       	MAX_LINE_LENGTH);
			strcpy(grid->test_data_file, test_filename);
			grid->traintype = TT;
		} else if (str_startswith(buffer, "p:")) {
			nr = all_doubles_str(buffer, 2, params);
			grid->ps = Calloc(double, nr);
			for (i=0; i<nr; i++)
				grid->ps[i] = params[i];
			grid->Np = nr;
		} else if (str_startswith(buffer, "lambda:")) {
			nr = all_doubles_str(buffer, 7, params);
			grid->lambdas = Calloc(double, nr);
			for (i=0; i<nr; i++)
				grid->lambdas[i] = params[i];
			grid->Nl = nr;
		} else if (str_startswith(buffer, "kappa:")) {
			nr = all_doubles_str(buffer, 6, params);
			grid->kappas = Calloc(double, nr);
			for (i=0; i<nr; i++)
				grid->kappas[i] = params[i];
			grid->Nk = nr;
		} else if (str_startswith(buffer, "epsilon:")) {
			nr = all_doubles_str(buffer, 8, params);
			grid->epsilons = Calloc(double, nr);
			for (i=0; i<nr; i++)
				grid->epsilons[i] = params[i];
			grid->Ne = nr;
		} else if (str_startswith(buffer, "weight:")) {
			nr = all_longs_str(buffer, 7, lparams);
			grid->weight_idxs = Calloc(int, nr);
			for (i=0; i<nr; i++)
				grid->weight_idxs[i] = lparams[i];
			grid->Nw = nr;
		} else if (str_startswith(buffer, "folds:")) {
			nr = all_longs_str(buffer, 6, lparams);
			grid->folds = lparams[0];
			if (nr > 1)
				fprintf(stderr, "Field \"folds\" only takes "
						"one value. Additional "
						"fields are ignored.\n");
		} else if (str_startswith(buffer, "repeats:")) {
			nr = all_longs_str(buffer, 8, lparams);
			grid->repeats = lparams[0];
			if (nr > 1)
				fprintf(stderr, "Field \"repeats\" only "
						"takes one value. Additional "
						"fields are ignored.\n");
		} else if (str_startswith(buffer, "kernel:")) {
			grid->kerneltype = parse_kernel_str(buffer);
		} else if (str_startswith(buffer, "gamma:")) {
			nr = all_doubles_str(buffer, 6, params);
			if (grid->kerneltype == K_LINEAR) {
				fprintf(stderr, "Field \"gamma\" ignored, "
						"linear kernel is used.\n");
				grid->Ng = 0;
				break;
			}
			grid->gammas = Calloc(double, nr);
			for (i=0; i<nr; i++)
				grid->gammas[i] = params[i];
			grid->Ng = nr;
		} else if (str_startswith(buffer, "coef:")) {
			nr = all_doubles_str(buffer, 5, params);
			if (grid->kerneltype == K_LINEAR ||
				grid->kerneltype == K_RBF) {
				fprintf(stderr, "Field \"coef\" ignored with "
						"specified kernel.\n");
				grid->Nc = 0;
				break;
			}
			grid->coefs = Calloc(double, nr);
			for (i=0; i<nr; i++)
				grid->coefs[i] = params[i];
			grid->Nc = nr;
		} else if (str_startswith(buffer, "degree:")) {
			nr = all_doubles_str(buffer, 7, params);
			if (grid->kerneltype != K_POLY) {
				fprintf(stderr, "Field \"degree\" ignored "
						"with specified kernel.\n");
				grid->Nd = 0;
				break;
			}
			grid->degrees = Calloc(double, nr);
			for (i=0; i<nr; i++)
				grid->degrees[i] = params[i];
			grid->Nd = nr;
		} else {
			fprintf(stderr, "Cannot find any parameters on line: "
					"%s\n", buffer);
		}
	}

	free(params);
	free(lparams);
	fclose(fid);
}
