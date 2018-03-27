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

#include "gensvm_checks.h"
#include "gensvm_cmdarg.h"
#include "gensvm_io.h"
#include "gensvm_train.h"
#include "gensvm_predict.h"

/**
 * Minimal number of command line arguments
 */
#define MINARGS 2

extern FILE *GENSVM_OUTPUT_FILE;
extern FILE *GENSVM_ERROR_FILE;

// function declarations
void exit_with_help(char **argv);
void parse_command_line(int argc, char **argv, struct GenModel *model,
		char **model_inputfile, char **training_inputfile,
		char **testing_inputfile, char **model_outputfile,
		char **prediction_outputfile);

/**
 * @brief Help function
 *
 * @details
 * Print help for this program and exit. Note that the VERSION is defined in 
 * the Makefile.
 *
 * @param[in] 	argv 	command line arguments
 *
 */
void exit_with_help(char **argv)
{
	printf("This is GenSVM, version %s.\n", VERSION_STRING);
	printf("Copyright (C) 2016, G.J.J. van den Burg.\n");
	printf("This program is free software, see the LICENSE file "
			"for details.\n\n");
	printf("Usage: %s [options] training_data [test_data]\n\n", argv[0]);
	printf("Options:\n");
	printf("--------\n");
	printf("-c coef              : coefficient for the polynomial and "
			"sigmoid kernel\n");
	printf("-d degree            : degree for the polynomial kernel\n");
	printf("-e epsilon           : set the value of the stopping "
			"criterion (epsilon > 0)\n");
	printf("-g gamma             : parameter for the rbf, polynomial or "
			"sigmoid kernel\n");
	printf("-h | -help           : print this help.\n");
	printf("-i max_iter          : maximum number of iterations to do.\n");
	printf("-k kappa             : set the value of kappa used in the "
			"Huber hinge (kappa > -1.0)\n");
	printf("-l lambda            : set the value of lambda "
			"(lambda > 0)\n");
	printf("-m model_output_file : write model output to file "
			"(not saved if no file provided)\n");
	printf("-o prediction_output : write predictions of test data to "
			"file (uses stdout if not provided)\n");
	printf("-p p-value           : set the value of p in the lp norm "
			"(1.0 <= p <= 2.0)\n");
	printf("-q                   : quiet mode (no output, not even "
			"errors!)\n");
	printf("-r rho               : choose the weigth specification "
			"(1 = unit, 2 = group)\n");
	printf("-s seed_model_file   : use previous model as seed for V\n");
	printf("-t type              : kerneltype (0=LINEAR, 1=POLY, 2=RBF, "
			"3=SIGMOID)\n");
	printf("-x                   : data files are in LibSVM/SVMlight "
			"format\n");
	printf("-z seed              : seed for the random number generator\n");
	printf("\n");

	exit(EXIT_FAILURE);
}

/**
 * @brief Main interface function for GenSVMtraintest
 *
 * @details
 * Main interface for the GenSVMtraintest commandline program.
 *
 * @param[in] 	argc 	number of command line arguments
 * @param[in] 	argv 	array of command line arguments
 *
 * @return 		exit status
 */
int main(int argc, char **argv)
{
	bool libsvm_format = false;
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
		exit_with_help(argv);

	// parse command line arguments
	parse_command_line(argc, argv, model, &model_inputfile,
		       	&training_inputfile, &testing_inputfile,
		       	&model_outputfile, &prediction_outputfile);
	libsvm_format = gensvm_check_argv(argc, argv, "-x");

	// read data from file
	if (libsvm_format)
		gensvm_read_data_libsvm(traindata, training_inputfile);
	else
		gensvm_read_data(traindata, training_inputfile);

	// check labels for consistency
	if (!gensvm_check_outcome_contiguous(traindata)) {
		err("[GenSVM Error]: Class labels should start from 1 and "
				"have no gaps. Please reformat your data.\n");
		exit(EXIT_FAILURE);
	}

	// save data filename to model
	model->data_file = Calloc(char, GENSVM_MAX_LINE_LENGTH);
	strcpy(model->data_file, training_inputfile);

	// check if we are sparse and want nonlinearity
	if (traindata->Z == NULL && model->kerneltype != K_LINEAR) {
		err("[GenSVM Warning]: Sparse matrices with nonlinear kernels "
				"are not yet supported. Dense matrices will "
				"be used.\n");
		traindata->RAW = gensvm_sparse_to_dense(traindata->spZ);
		traindata->Z = traindata->RAW;
		gensvm_free_sparse(traindata->spZ);
	}

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
		// read the test data
		if (libsvm_format)
			gensvm_read_data_libsvm(testdata, testing_inputfile);
		else
			gensvm_read_data(testdata, testing_inputfile);

		// check if we are sparse and want nonlinearity
		if (testdata->Z == NULL && model->kerneltype != K_LINEAR) {
			err("[GenSVM Warning]: Sparse matrices with nonlinear "
					"kernels are not yet supported. Dense "
					"matrices will be used.\n");
			testdata->Z = gensvm_sparse_to_dense(testdata->spZ);
			gensvm_free_sparse(testdata->spZ);
		}

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

/**
 * @brief Exit with warning about invalid parameter value.
 *
 * @param[in] 	label 	name of the parameter
 * @param[in] 	argv 	command line arguments
 */
void exit_invalid_param(const char *label, char **argv)
{
	fprintf(stderr, "Invalid parameter value for %s.\n\n", label);
	exit_with_help(argv);
}

/**
 * @brief Parse the command line arguments
 *
 * @details
 * For a full overview of the command line arguments and their meaning see 
 * exit_with_help(). This function furthermore sets the default output streams 
 * to stdout/stderr, and initializes the kernel parameters if none are 
 * supplied: gamma = 1.0, degree = 2.0, coef = 0.0.
 *
 * @param[in] 	 argc 			number of command line arguments
 * @param[in] 	 argv 			array of command line arguments
 * @param[in] 	 model 			initialized GenModel struct
 * @param[out] 	 model_inputfile 	filename for the seed model
 * @param[out] 	 training_inputfile 	filename for the training data
 * @param[out] 	 testing_inputfile 	filename for the test data
 * @param[out] 	 model_outputfile 	filename for the output model
 * @param[out] 	 prediction_outputfile 	filename for the predictions
 *
 */
void parse_command_line(int argc, char **argv, struct GenModel *model,
		char **model_inputfile, char **training_inputfile,
	       	char **testing_inputfile, char **model_outputfile,
	       	char **prediction_outputfile)
{
	int i;

	GENSVM_OUTPUT_FILE = stdout;
	GENSVM_ERROR_FILE = stderr;

	// parse options
	// note: flags that don't have an argument should decrement i
	for (i=1; i<argc; i++) {
		if (argv[i][0] != '-') break;
		if (++i>=argc) {
			exit_with_help(argv);
		}
		switch (argv[i-1][1]) {
			case 'c':
				model->coef = atof(argv[i]);
				break;
			case 'd':
				model->degree = atof(argv[i]);
				break;
			case 'e':
				model->epsilon = atof(argv[i]);
				if (model->epsilon <= 0)
					exit_invalid_param("epsilon", argv);
				break;
			case 'g':
				model->gamma = atof(argv[i]);
				break;
			case 'i':
				model->max_iter = atoi(argv[i]);
				break;
			case 'k':
				model->kappa = atof(argv[i]);
				if (model->kappa <= -1.0)
					exit_invalid_param("kappa", argv);
				break;
			case 'l':
				model->lambda = atof(argv[i]);
				if (model->lambda <= 0)
					exit_invalid_param("lambda", argv);
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
				if (model->p < 1.0 || model->p > 2.0)
					exit_invalid_param("p", argv);
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
			case 'x':
				i--;
				break;
			case 'z':
				model->seed = atoi(argv[i]);
				break;
			default:
				// this one should always print explicitly to 
				// stderr, even if '-q' is supplied, because 
				// otherwise you can't debug cmdline flags.
				fprintf(stderr, "Unknown option: -%c\n",
						argv[i-1][1]);
				exit_with_help(argv);
		}
	}
	if (i >= argc)
		exit_with_help(argv);

	(*training_inputfile) = Malloc(char, strlen(argv[i])+1);
	strcpy((*training_inputfile), argv[i]);
	if (i+2 == argc) {
		(*testing_inputfile) = Malloc(char, strlen(argv[i])+1);
		strcpy((*testing_inputfile), argv[i+1]);
	}
}
