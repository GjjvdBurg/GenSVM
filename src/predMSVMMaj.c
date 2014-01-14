/**
 * @file predMSVMMaj.c
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Command line interface for predicting class labels
 *
 * @details
 * This is a command line program for predicting the class labels or
 * determining the predictive performance of a pre-determined model on a given
 * test dataset. The predictive performance can be written to the screen or
 * the predicted class labels can be written to a specified output file. This
 * is done using msvmmaj_write_predictions().
 *
 * The specified model file must follow the specification given in
 * msvmmaj_write_model().
 *
 * For usage information, see the program help function.
 *
 */

#include "msvmmaj.h"
#include "msvmmaj_init.h"
#include "msvmmaj_pred.h"
#include "util.h"

#define MINARGS 3

extern FILE *MSVMMAJ_OUTPUT_FILE;

// function declarations
void exit_with_help();
void parse_command_line(int argc, char **argv, 
		char *input_filename, char *output_filename, 
		char *model_filename);

/**
 * @brief Help function
 */
void exit_with_help()
{
	printf("This is MSVMMaj, version %1.1f\n\n", VERSION);
	printf("Usage: predMSVMMaj [options] test_data_file model_file\n");
	printf("Options:\n");
	printf("-o output_file : write output to file\n");
	printf("-q : quiet mode (no output)\n");
	exit(0);
}

/**
 * @brief Main interface function for predMSVMMaj
 *
 * @details
 * Main interface for the command line program. A given model file is read and
 * a test dataset is initialized from the given data. The predictive
 * performance (hitrate) of the model on the test set is printed to the output
 * stream (default = stdout). If an output file is specified the predictions
 * are written to the file.
 *
 * @todo
 * Ensure that the program can read model files without class labels
 * specified. In that case no prediction accuracy is printed to the screen.
 *
 * @param[in] 	argc 	number of command line arguments
 * @param[in] 	argv 	array of command line arguments
 *
 */
int main(int argc, char **argv)
{
	long *predy;
	double performance;

	char input_filename[MAX_LINE_LENGTH];
	char model_filename[MAX_LINE_LENGTH];
	char output_filename[MAX_LINE_LENGTH];;

	if (argc < MINARGS || msvmmaj_check_argv(argc, argv, "-help") 
			|| msvmmaj_check_argv_eq(argc, argv, "-h") )
		exit_with_help();
	parse_command_line(argc, argv, input_filename, output_filename,
			model_filename);

	// read the data and model
	struct MajModel *model = msvmmaj_init_model();
	struct MajData *data = msvmmaj_init_data();
	msvmmaj_read_data(data, input_filename);
	msvmmaj_read_model(model, model_filename);

	// check if the number of attributes in data equals that in model
	if (data->m != model->m) {
		fprintf(stderr, "Error: number of attributes in data (%li) "
				"does not equal the number of attributes in "
				"model (%li)\n", data->m, model->m);
		exit(1);
	} else if (data->K != model->K) {
		fprintf(stderr, "Error: number of classes in data (%li) "
				"does not equal the number of classes in "
				"model (%li)\n", data->K, model->K);
		exit(1);
	}

	// predict labels and performance if test data has labels
	predy = Calloc(long, data->n);
	msvmmaj_predict_labels(data, model, predy);
	if (data->y != NULL) {
		performance = msvmmaj_prediction_perf(data, predy);
		note("Predictive performance: %3.2f%%\n", performance);
	}

	// if output file is specified, write predictions to it
	if (msvmmaj_check_argv_eq(argc, argv, "-o")) {
		msvmmaj_write_predictions(data, predy, output_filename);
		note("Predictions written to: %s\n", output_filename);
	}

	// free the model, data, and predictions
	msvmmaj_free_model(model);
	msvmmaj_free_data(data);
	free(predy);

	return 0;
}

/**
 * @brief Parse command line arguments
 *
 * @details
 * Read the data filename and model filename from the command line arguments.
 * If specified, also read the output filename. If the quiet flag is given,
 * set the global output stream to NULL. On error, exit_with_help().
 *
 * @param[in] 	argc 			number of command line arguments
 * @param[in] 	argv 			array of command line arguments
 * @param[in] 	input_filename 		pre-allocated array for the input 
 * 					filename
 * @param[in] 	output_filename 	pre-allocated array for the output
 * 					filename
 * @param[in] 	model_filename 		pre-allocated array for the model
 * 					filename
 *
 */
void parse_command_line(int argc, char **argv, char *input_filename,
		char *output_filename, char *model_filename)
{
	int i;

	MSVMMAJ_OUTPUT_FILE = stdout;

	for (i=1; i<argc; i++) {
		if (argv[i][0] != '-') break;
		if (++i >= argc)
			exit_with_help();
		switch (argv[i-1][1]) {
			case 'o':
				strcpy(output_filename, argv[i]);
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

	if (i >= argc)
		exit_with_help();

	strcpy(input_filename, argv[i]);
	i++;
	strcpy(model_filename, argv[i]);
}
