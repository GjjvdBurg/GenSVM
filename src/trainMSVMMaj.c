#include <time.h>
#include <math.h>

#include "libMSVMMaj.h"
#include "msvmmaj_train.h"
#include "util.h"
#include "MSVMMaj.h"

#define MINARGS 2

extern FILE *MSVMMAJ_OUTPUT_FILE;

void print_null(const char *s) {}
void exit_with_help();
void parse_command_line(int argc, char **argv, struct MajModel *model, 
		char *input_filename, char *output_filename, char *model_filename);

void exit_with_help()
{
	printf("This is MSVMMaj, version %1.1f\n\n", VERSION);
	printf("Usage: trainMSVMMaj [options] training_data_file\n");
	printf("Options:\n");
	printf("-e epsilon : set the value of the stopping criterion\n");
	printf("-h | -help : print this help.\n");
	printf("-k kappa : set the value of kappa used in the Huber hinge\n");
	printf("-l lambda : set the value of lambda (lambda > 0)\n");
	printf("-m model_file : use previous model as seed for W and t\n");
	printf("-o output_file : write output to file\n");
	printf("-p p-value : set the value of p in the lp norm (1.0 <= p <= 2.0)\n");
	printf("-q : quiet mode (no output)\n");
	printf("-r rho : choose the weigth specification (1 = unit, 2 = group)\n");

	exit(0);
}

/*
	Main
*/
int main(int argc, char **argv)
{
	char input_filename[MAX_LINE_LENGTH];
	char model_filename[MAX_LINE_LENGTH];
	char output_filename[MAX_LINE_LENGTH];

	struct MajModel *model = Malloc(struct MajModel, 1);
	struct MajData *data = Malloc(struct MajData, 1);

	if (argc < MINARGS || msvmmaj_check_argv(argc, argv, "-help") 
			|| msvmmaj_check_argv_eq(argc, argv, "-h") ) 
		exit_with_help();
	parse_command_line(argc, argv, model, input_filename, output_filename, model_filename);

	// read data file
	msvmmaj_read_data(data, input_filename);

	// copy dataset parameters to model	
	model->n = data->n;
	model->m = data->m;
	model->K = data->K;
	model->data_file = input_filename;

	// allocate model and initialize weights
	msvmmaj_allocate_model(model);
	msvmmaj_initialize_weights(data, model);

	srand(time(NULL));

	if (msvmmaj_check_argv_eq(argc, argv, "-m")) {
		struct MajModel *seed_model = Malloc(struct MajModel, 1);
		msvmmaj_read_model(seed_model, model_filename);
		msvmmaj_seed_model_V(seed_model, model);
		msvmmaj_free_model(seed_model);
	} else {
		msvmmaj_seed_model_V(NULL, model);
	}
	// initialize kernel (if necessary)
	// msvmmaj_make_kernel(model, data);

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

void parse_command_line(int argc, char **argv, struct MajModel *model, 
		char *input_filename, char *output_filename, char *model_filename)
{
	int i;

	// default values
	model->p = 1.0;
	model->lambda = pow(2, -8.0);
	model->epsilon = 1e-6;
	model->kappa = 0.0;
	model->weight_idx = 1;
	
	MSVMMAJ_OUTPUT_FILE = stdout;

	// parse options
	for (i=1; i<argc; i++) {
		if (argv[i][0] != '-') break;
		if (++i>=argc) {
			exit_with_help();
		}
		switch (argv[i-1][1]) {
			case 'e':
				model->epsilon = atof(argv[i]);
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
			case 'q':
				MSVMMAJ_OUTPUT_FILE = NULL;
				i--;
				break;
			default:
				fprintf(stderr, "Unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
		}
	}

	// read input filename
	if (i >= argc)
		exit_with_help();

	strcpy(input_filename, argv[i]);
}

