/**
 * @file GenSVMtraintest.c
 * @author G.J.J. van den Burg
 * @date 2015-02-01
 * @brief Command line interface for training and testing with a GenSVM model
 *
 * @details
 * This is a command line program for training and testing on a single model
 * with specified model parameters.
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

#include "gensvm_cmdarg.h"
#include "gensvm_io.h"
#include "gensvm_train.h"
#include "gensvm_predict.h"

#define MINARGS 2

extern FILE *GENSVM_OUTPUT_FILE;
extern FILE *GENSVM_ERROR_FILE;

// function declarations
void exit_with_help();
void parse_command_line(int argc, char **argv, struct GenModel *model,
		char **model_inputfile, char **training_inputfile,
		char **testing_inputfile, char **model_outputfile,
		char **prediction_outputfile);

void exit_with_help()
{
	printf("This is GenSVM, version %1.1f\n\n", VERSION);
	printf("Usage: ./gensvm [options] training_data [test_data]\n");
	printf("Options:\n");
	printf("-c coef : coefficient for the polynomial and sigmoid kernel\n");
	printf("-d degree : degree for the polynomial kernel\n");
	printf("-e epsilon : set the value of the stopping criterion\n");
	printf("-g gamma : parameter for the rbf, polynomial or sigmoid "
			"kernel\n");
	printf("-h | -help : print this help.\n");
	printf("-k kappa : set the value of kappa used in the Huber hinge\n");
	printf("-l lambda : set the value of lambda (lambda > 0)\n");
	printf("-s seed_model_file : use previous model as seed for V\n");
	printf("-m model_output_file : write model output to file\n");
	printf("-o prediction_output : write predictions of test data to "
			"file\n");
	printf("-p p-value : set the value of p in the lp norm "
			"(1.0 <= p <= 2.0)\n");
	printf("-q : quiet mode (no output, not even errors!)\n");
	printf("-r rho : choose the weigth specification (1 = unit, 2 = "
			"group)\n");
	printf("-t type: kerneltype (0=LINEAR, 1=POLY, 2=RBF, 3=SIGMOID)\n");

	exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
	long i, *predy = NULL;
	double performance;

	char *training_inputfile = NULL,
	     *testing_inputfile = NULL,
	     *model_inputfile = NULL,
	     *model_outputfile = NULL,
	     *prediction_outputfile = NULL;

	struct GenModel *model = gensvm_init_model();
	struct GenModel *seed_model = NULL;
	struct GenData *traindata = gensvm_init_data();
	struct GenData *testdata = gensvm_init_data();

	if (argc < MINARGS || gensvm_check_argv(argc, argv, "-help")
		       	|| gensvm_check_argv_eq(argc, argv, "-h"))
		exit_with_help();
	parse_command_line(argc, argv, model, &model_inputfile,
		       	&training_inputfile, &testing_inputfile,
		       	&model_outputfile, &prediction_outputfile);

	// read data from files
	gensvm_read_data(traindata, training_inputfile);
	model->data_file = Calloc(char, MAX_LINE_LENGTH);
	strcpy(model->data_file, training_inputfile);

	// seed the random number generator
	srand(time(NULL));

	// load a seed model from file if it is specified
	if (gensvm_check_argv_eq(argc, argv, "-s")) {
		seed_model = gensvm_init_model();
		gensvm_read_model(seed_model, model_inputfile);
	}

	// train the GenSVM model
	gensvm_train(model, traindata, seed_model);

	// if we also have a test set, predict labels and write to predictions
	// to an output file if specified
	if (testing_inputfile != NULL) {
		gensvm_read_data(testdata, testing_inputfile);
		gensvm_kernel_postprocess(model, traindata, testdata);

		// predict labels
		predy = Calloc(long, testdata->n);
		gensvm_predict_labels(testdata, model, predy);

		if (testdata->y != NULL) {
			performance = gensvm_prediction_perf(testdata, predy);
			note("Predictive performance: %3.2f%%\n", performance);
		}

		// if output file is specified, write predictions to it
		if (gensvm_check_argv_eq(argc, argv, "-o")) {
			gensvm_write_predictions(testdata, predy,
				       	prediction_outputfile);
			note("Prediction written to: %s\n",
				       	prediction_outputfile);
		} else {
			for (i=0; i<testdata->n; i++)
				printf("%li ", predy[i]);
			printf("\n");
		}
	}

	// write model to output file if necessary
	if (gensvm_check_argv_eq(argc, argv, "-m")) {
		gensvm_write_model(model, model_outputfile);
		note("Model written to: %s\n", model_outputfile);
	}

	// free everything
	gensvm_free_model(model);
	gensvm_free_model(seed_model);
	gensvm_free_data(traindata);
	gensvm_free_data(testdata);
	free(training_inputfile);
	free(testing_inputfile);
	free(model_inputfile);
	free(model_outputfile);
	free(prediction_outputfile);

	free(predy);

	return 0;
}

void parse_command_line(int argc, char **argv, struct GenModel *model,
		char **model_inputfile, char **training_inputfile,
	       	char **testing_inputfile, char **model_outputfile,
	       	char **prediction_outputfile)
{
	int i;
	double gamma = 1.0,
	       degree = 2.0,
	       coef = 0.0;

	GENSVM_OUTPUT_FILE = stdout;
	GENSVM_ERROR_FILE = stderr;

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
			case 's':
				(*model_inputfile) = Malloc(char,
					       	strlen(argv[i])+1);
				strcpy((*model_inputfile), argv[i]);
				break;
			case 'm':
				(*model_outputfile) = Malloc(char,
						strlen(argv[i])+1);
				strcpy((*model_outputfile), argv[i]);
				break;
			case 'o':
				(*prediction_outputfile) = Malloc(char,
						strlen(argv[i])+1);
				strcpy((*prediction_outputfile), argv[i]);
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
				GENSVM_OUTPUT_FILE = NULL;
				GENSVM_ERROR_FILE = NULL;
				i--;
				break;
			default:
				// this one should always print explicitly to 
				// stderr, even if '-q' is supplied, because 
				// otherwise you can't debug cmdline flags.
				fprintf(stderr, "Unknown option: -%c\n",
						argv[i-1][1]);
				exit_with_help();
		}
	}
	if (i >= argc)
		exit_with_help();

	(*training_inputfile) = Malloc(char, strlen(argv[i])+1);
	strcpy((*training_inputfile), argv[i]);
	if (i+2 == argc) {
		(*testing_inputfile) = Malloc(char, strlen(argv[i])+1);
		strcpy((*testing_inputfile), argv[i+1]);
	}

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
