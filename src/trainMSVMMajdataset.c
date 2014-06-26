/**
 * @file trainMSVMMajdataset.c
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Command line interface for the grid search program
 *
 * @details
 * This is a command line interface to the parameter grid search functionality
 * of the algorithm. The grid search is specified in a separate file, thereby
 * reducing the number of command line arguments. See
 * read_training_from_file() for documentation on the training file. 
 *
 * The program runs a grid search as specified in the training file. If
 * desired the grid search can incorporate consistency checks to find the
 * configuration among the best configurations which scores consistently high.
 * All output is written to stdout, unless the quiet mode is specified.
 *
 * For further usage information, see the program help function.
 *
 */

#include <time.h>

#include "crossval.h"
#include "msvmmaj.h"
#include "msvmmaj_io.h"
#include "msvmmaj_init.h"
#include "msvmmaj_pred.h"
#include "msvmmaj_train.h"
#include "msvmmaj_train_dataset.h"
#include "strutil.h"
#include "util.h"

#define MINARGS 2

extern FILE *MSVMMAJ_OUTPUT_FILE;

// function declarations
void exit_with_help();
void parse_command_line(int argc, char **argv, char *input_filename);
void read_training_from_file(char *input_filename, struct Training *training);

/**
 * @brief Help function
 */
void exit_with_help()
{
	printf("This is MSVMMaj, version %1.1f\n\n", VERSION);
	printf("Usage: trainMSVMMajdataset [options] training_file\n");
	printf("Options:\n");
	printf("-h | -help : print this help.\n");
	printf("-q : quiet mode (no output)\n");

	exit(0);
}

/**
 * @brief Main interface function for trainMSVMMajdataset
 *
 * @details
 * Main interface for the command line program. A given training file which
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
	
	struct Training *training = Malloc(struct Training, 1);	
	struct MajData *train_data = Malloc(struct MajData, 1);
	struct MajData *test_data = Malloc(struct MajData, 1);

	if (argc < MINARGS || msvmmaj_check_argv(argc, argv, "-help") 
			|| msvmmaj_check_argv_eq(argc, argv, "-h") )
		exit_with_help();
	parse_command_line(argc, argv, input_filename);

	training->repeats = 0;
	note("Reading training file\n");	
	read_training_from_file(input_filename, training);

	note("Reading data from %s\n", training->train_data_file);
	msvmmaj_read_data(train_data, training->train_data_file);
	if (training->traintype == TT) {
		note("Reading data from %s\n", training->test_data_file);
		msvmmaj_read_data(test_data, training->test_data_file);
	}

	note("Creating queue\n");
	struct Queue *q = Malloc(struct Queue, 1);
	make_queue(training, q, train_data, test_data);

	srand(time(NULL));

	note("Starting training\n");
	if (training->traintype == TT)
		start_training_tt(q);
	else
		start_training_cv(q);
	note("Training finished\n");

	if (training->repeats > 0) {
		consistency_repeats(q, training->repeats, training->traintype);
	}

	printf("here 0\n");
	free_queue(q);
	printf("here 1\n");
	free(training);
	printf("here 2\n");
	msvmmaj_free_data(train_data);
	printf("here 3\n");
	msvmmaj_free_data(test_data);
	printf("here 4\n");

	note("Done.\n");
	return 0;
}

/**
 * @brief Parse command line arguments
 *
 * @details
 * Few arguments can be supplied to the command line. Only quiet mode can be
 * specified, or help can be requested. The filename of the training file is
 * read from the arguments. Parsing of the training file is done separately in
 * read_training_from_file().
 *
 * @param[in] 	argc 		number of command line arguments
 * @param[in] 	argv 		array of command line arguments
 * @param[in] 	input_filename 	pre-allocated buffer for the training
 * 				filename.
 *
 */
void parse_command_line(int argc, char **argv, char *input_filename)
{
	int i;

	MSVMMAJ_OUTPUT_FILE = stdout;

	for (i=1; i<argc; i++) {
		if (argv[i][0] != '-') break;
		if (++i>=argc)
			exit_with_help();
		switch (argv[i-1][1]) {
			case 'q':
				MSVMMAJ_OUTPUT_FILE = NULL;
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
 * @brief Read the Training struct from file
 *
 * @details
 * Read the Training struct from a file. The training file follows a specific
 * format specified in @ref spec_training_file. 
 *
 * Commonly used string functions in this function are all_doubles_str() and
 * all_longs_str(). 
 *
 * @param[in] 	input_filename 	filename of the training file
 * @param[in] 	training 	Training structure to place the parsed
 * 				parameter grid.
 *
 */
void read_training_from_file(char *input_filename, struct Training *training)
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
		fprintf(stderr, "Error opening training file %s\n", 
				input_filename);
		exit(1);
	}
	training->traintype = CV;
	while ( fgets(buffer, MAX_LINE_LENGTH, fid) != NULL ) {
		Memset(params, double,  MAX_LINE_LENGTH);
		Memset(lparams, long, MAX_LINE_LENGTH);
		if (str_startswith(buffer, "train:")) {
			sscanf(buffer, "train: %s\n", train_filename);
			training->train_data_file = Calloc(char, 
					MAX_LINE_LENGTH);
			strcpy(training->train_data_file, train_filename);
		} else if (str_startswith(buffer, "test:")) {
			sscanf(buffer, "test: %s\n", test_filename);
			training->test_data_file = Calloc(char,
				       	MAX_LINE_LENGTH);
			strcpy(training->test_data_file, test_filename);
			training->traintype = TT;
		} else if (str_startswith(buffer, "p:")) {
			nr = all_doubles_str(buffer, 2, params);
			training->ps = Calloc(double, nr);
			for (i=0; i<nr; i++)
				training->ps[i] = params[i];
			training->Np = nr;
		} else if (str_startswith(buffer, "lambda:")) {
			nr = all_doubles_str(buffer, 7, params);
			training->lambdas = Calloc(double, nr);
			for (i=0; i<nr; i++)
				training->lambdas[i] = params[i];
			training->Nl = nr;
		} else if (str_startswith(buffer, "kappa:")) {
			nr = all_doubles_str(buffer, 6, params);
			training->kappas = Calloc(double, nr);
			for (i=0; i<nr; i++)
				training->kappas[i] = params[i];
			training->Nk = nr;
		} else if (str_startswith(buffer, "epsilon:")) {
			nr = all_doubles_str(buffer, 8, params);
			training->epsilons = Calloc(double, nr);
			for (i=0; i<nr; i++)
				training->epsilons[i] = params[i];
			training->Ne = nr;
		} else if (str_startswith(buffer, "weight:")) {
			nr = all_longs_str(buffer, 7, lparams);
			training->weight_idxs = Calloc(int, nr);
			for (i=0; i<nr; i++)
				training->weight_idxs[i] = lparams[i];
			training->Nw = nr;
		} else if (str_startswith(buffer, "folds:")) {
			nr = all_longs_str(buffer, 6, lparams);
			training->folds = lparams[0];
			if (nr > 1)
				fprintf(stderr, "Field \"folds\" only takes "
						"one value. Additional "
						"fields are ignored.\n");
		} else if (str_startswith(buffer, "repeats:")) {
			nr = all_longs_str(buffer, 8, lparams);
			training->repeats = lparams[0];
			if (nr > 1)
				fprintf(stderr, "Field \"repeats\" only "
						"takes one value. Additional "
						"fields are ignored.\n");
		} else if (str_startswith(buffer, "kernel:")) {
			training->kerneltype = parse_kernel_str(buffer);
		} else if (str_startswith(buffer, "gamma:")) {
			nr = all_doubles_str(buffer, 6, params);
			if (training->kerneltype == K_LINEAR) {
				fprintf(stderr, "Field \"gamma\" ignored, "
						"linear kernel is used.\n");
				training->Ng = 0;
				break;
			}
			training->gammas = Calloc(double, nr);
			for (i=0; i<nr; i++)
				training->gammas[i] = params[i];
			training->Ng = nr;
		} else if (str_startswith(buffer, "coef:")) {
			nr = all_doubles_str(buffer, 5, params);
			if (training->kerneltype == K_LINEAR ||
				training->kerneltype == K_RBF) {
				fprintf(stderr, "Field \"coef\" ignored with "
						"specified kernel.\n");
				training->Nc = 0;
				break;
			}
			training->coefs = Calloc(double, nr);
			for (i=0; i<nr; i++)
				training->coefs[i] = params[i];
			training->Nc = nr;
		} else if (str_startswith(buffer, "degree:")) {
			nr = all_doubles_str(buffer, 7, params);
			if (training->kerneltype != K_POLY) {
				fprintf(stderr, "Field \"degree\" ignored "
						"with specified kernel.\n");
				training->Nd = 0;
				break;
			}
			training->degrees = Calloc(double, nr);
			for (i=0; i<nr; i++)
				training->degrees[i] = params[i];
			training->Nd = nr;
		} else {
			fprintf(stderr, "Cannot find any parameters on line: "
					"%s\n", buffer);
		}
	}

	free(params);
	free(lparams);
	fclose(fid);
}
