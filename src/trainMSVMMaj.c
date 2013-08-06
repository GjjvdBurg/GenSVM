#include "libMSVMMaj.h"

#define MINARGS 2

void print_null(const char *s) {}
void exit_with_help();
void parse_command_line(int argc, char **argv, struct Model *model, 
		char *input_filename, char *output_filename);

void exit_with_help()
{
	printf("This is MSVMMaj, version %1.1f\n\n", VERSION);
	printf("Usage: trainMSVMMaj [options] training_data_file\n");
	printf("Options:\n");
	printf("-c folds : perform cross validation with given number of folds\n");
	printf("-e epsilon : set the value of the stopping criterion\n");
	printf("-h | -help : print this help.\n");
	printf("-k kappa : set the value of kappa used in the Huber hinge\n");
	printf("-l lambda : set the value of lambda (lambda > 0)\n");
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

	struct Model *model = Malloc(struct Model, 1);
	struct Data *data = Malloc(struct Data, 1);

	if (argc < MINARGS || check_argv(argc, argv, "-help") || check_argv_eq(argc, argv, "-h") ) 
		exit_with_help();
	parse_command_line(argc, argv, model, input_filename, model_filename);

	// read data file
	read_data(data, input_filename);

	// copy dataset parameters to model	
	model->n = data->n;
	model->m = data->m;
	model->K = data->K;
	model->data_file = input_filename;

	// allocate model and initialize weights
	allocate_model(model);
	initialize_weights(data, model);

	// start training
	main_loop(model, data);

	// write_model to file
	if (check_argv_eq(argc, argv, "-o")) {
		write_model(model, model_filename);
		info("Output written to %s\n", model_filename);
	}

	// free model and data
	free_model(model);
	free_data(data);
	
	return 0;
}

void parse_command_line(int argc, char **argv, struct Model *model, 
		char *input_filename, char *output_filename)
{
	int i;
	void (*print_func)(const char*) = NULL;

	// default values
	model->p = 1.0;
	model->lambda = pow(2, -8.0);
	model->epsilon = 1e-6;
	model->kappa = 0.0;
	model->weight_idx = 1;

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
				print_func = &print_null;
				i--;
				break;
			default:
				fprintf(stderr, "Unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
		}
	}

	// set print function
	set_print_string_function(print_func);

	// read input filename
	if (i >= argc)
		exit_with_help();

	strcpy(input_filename, argv[i]);
}

