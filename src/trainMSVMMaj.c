/**
 * @file trainMSVMMaj.c
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Command line interface for training a single model with MSVMMaj
 *
 * @details
 * This is a command line program for training a single model on a given
 * dataset. To run a grid search over a number of parameter configurations, 
 * see trainMSVMMajdataset.c. 
 *
 */

#include <time.h>
#include <math.h>

#include "msvmmaj_kernel.h"
#include "libMSVMMaj.h"
#include "msvmmaj.h"
#include "msvmmaj_io.h"
#include "msvmmaj_init.h"
#include "msvmmaj_train.h"
#include "util.h"

#define MINARGS 2

extern FILE *MSVMMAJ_OUTPUT_FILE;

// function declarations
void exit_with_help();
void parse_command_line(int argc, char **argv, struct MajModel *model, 
		char *input_filename, char *output_filename, char *model_filename);

/**
 * @brief Help function
 */
void exit_with_help()
{
	printf("This is MSVMMaj, version %1.1f\n\n", VERSION);
	printf("Usage: trainMSVMMaj [options] training_data_file\n");
	printf("Options:\n");
	printf("-c coef : coefficient for the polynomial and sigmoid kernel\n");
	printf("-d degree : degree for the polynomial kernel\n");
	printf("-e epsilon : set the value of the stopping criterion\n");
	printf("-g gamma : parameter for the rbf, polynomial or sigmoid "
			"kernel\n");
	printf("-h | -help : print this help.\n");
	printf("-k kappa : set the value of kappa used in the Huber hinge\n");
	printf("-l lambda : set the value of lambda (lambda > 0)\n");
	printf("-m model_file : use previous model as seed for W and t\n");
	printf("-o output_file : write output to file\n");
	printf("-p p-value : set the value of p in the lp norm "
			"(1.0 <= p <= 2.0)\n");
	printf("-q : quiet mode (no output)\n");
	printf("-r rho : choose the weigth specification (1 = unit, 2 = "
			"group)\n");
	printf("-t type: kerneltype (LINEAR=0, POLY=1, RBF=2, SIGMOID=3)\n");

	exit(0);
}

/**
 * @brief Main interface function for trainMSVMMaj
 *
 * @details
 * Main interface for the command line program. A given dataset file is read
 * and a MSVMMaj model is trained on this data. By default the progress of the
 * computations are written to stdout. See for full options of the program the
 * help function.
 *
 * @param[in] 	argc 	number of command line arguments
 * @param[in] 	argv 	array of command line arguments
 *
 */
int main(int argc, char **argv)
{
	char input_filename[MAX_LINE_LENGTH];
	char model_filename[MAX_LINE_LENGTH];
	char output_filename[MAX_LINE_LENGTH];

	struct MajModel *model = msvmmaj_init_model();
	struct MajData *data = msvmmaj_init_data();

	if (argc < MINARGS || msvmmaj_check_argv(argc, argv, "-help") 
			|| msvmmaj_check_argv_eq(argc, argv, "-h") ) 
		exit_with_help();
	parse_command_line(argc, argv, model, input_filename, 
			output_filename, model_filename);

	// read data file
	msvmmaj_read_data(data, input_filename);

	// copy dataset parameters to model	
	model->n = data->n;
	model->m = data->m;
	model->K = data->K;
	model->data_file = input_filename;
	
	// allocate model
	msvmmaj_allocate_model(model);
	
	// initialize kernel (if necessary)
	msvmmaj_make_kernel(model, data);

	// reallocate model and initialize weights
	msvmmaj_reallocate_model(model, data->n, data->m);
	msvmmaj_initialize_weights(data, model);

	// seed the random number generator (only place in programs is in
	// command line interfaces)
	//srand(time(NULL));
	srand(123123);	

	if (msvmmaj_check_argv_eq(argc, argv, "-m")) {
		struct MajModel *seed_model = msvmmaj_init_model();
		msvmmaj_read_model(seed_model, model_filename);
		msvmmaj_seed_model_V(seed_model, model, data);
		msvmmaj_free_model(seed_model);
	} else {
		msvmmaj_seed_model_V(NULL, model, data);
	}

	// start training
	msvmmaj_optimize(model, data);

	// write_model to file
	if (msvmmaj_check_argv_eq(argc, argv, "-o")) {
		msvmmaj_write_model(model, output_filename);
		note("Output written to %s\n", output_filename);
	}

	// free model and data
	msvmmaj_free_model(model);
	msvmmaj_free_data(data);
	
	return 0;
}

/**
 * @brief Parse command line arguments
 *
 * @details
 * Process the command line arguments for the model parameters, and record
 * them in the specified MajModel. An input filename for the dataset is read
 * and if specified an output filename and a model filename for the seed
 * model.
 *
 * @param[in] 	argc 			number of command line arguments
 * @param[in] 	argv 			array of command line arguments
 * @param[in] 	model 			initialized model
 * @param[in] 	input_filename 		pre-allocated buffer for the input
 * 					filename
 * @param[in] 	output_filename 	pre-allocated buffer for the output
 * 					filename
 * @param[in] 	model_filename 		pre-allocated buffer for the model
 * 					filename
 *
 */
void parse_command_line(int argc, char **argv, struct MajModel *model, 
		char *input_filename, char *output_filename, char *model_filename)
{
	int i;
	double gamma = 1.0, 
	       degree = 2.0, 
	       coef = 0.0;

	MSVMMAJ_OUTPUT_FILE = stdout;

	// parse options
	for (i=1; i<argc; i++) {
		if (argv[i][0] != '-') break;
		if (++i>=argc) {
			exit_with_help();
		}
		switch (argv[i-1][1]) {
			case 'c':
				coef = atof(argv[i]);
				break;
			case 'd':
				degree = atof(argv[i]);
				break;
			case 'e':
				model->epsilon = atof(argv[i]);
				break;
			case 'g':
				gamma = atof(argv[i]);
				break;
			case 'k':
				model->kappa = atof(argv[i]);
				break;
			case 'l':
				model->lambda = atof(argv[i]);
				break;
			case 'm':
				strcpy(model_filename, argv[i]);
				break;
			case 'o':
				strcpy(output_filename, argv[i]);
				break;
			case 'p':
				model->p = atof(argv[i]);
				break;
			case 'r':
				model->weight_idx = atoi(argv[i]);
				break;
			case 't':
				model->kerneltype = atoi(argv[i]);
				break;
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
	
	// read input filename
	if (i >= argc)
		exit_with_help();

	strcpy(input_filename, argv[i]);

	// set kernel parameters
	switch (model->kerneltype) {
		case K_LINEAR:
			break;
		case K_POLY:
			model->kernelparam = Calloc(double, 3);
			model->kernelparam[0] = gamma;
			model->kernelparam[1] = coef;
			model->kernelparam[2] = degree;
			break;
		case K_RBF:
			model->kernelparam = Calloc(double, 1);
			model->kernelparam[0] = gamma;
			break;
		case K_SIGMOID:
			model->kernelparam = Calloc(double, 1);
			model->kernelparam[0] = gamma;
			model->kernelparam[1] = coef;
	}
}
