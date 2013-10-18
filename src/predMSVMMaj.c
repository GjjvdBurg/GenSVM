#include "msvmmaj_pred.h"
#include "MSVMMaj.h"
#include "util.h"

#define MINARGS 3

extern FILE *MSVMMAJ_OUTPUT_FILE;

void print_null(const char *s) {}
void exit_with_help();
void parse_command_line(int argc, char **argv, struct MajModel *model,
		char *input_filename, char *output_filename, 
		char *model_filename);

void exit_with_help()
{
	printf("This is MSVMMaj, version %1.1f\n\n", VERSION);
	printf("Usage: predMSVMMaj [options] test_data_file model_file\n");
	printf("Options:\n");
	printf("-o output_file : write output to file\n");
	printf("-q : quiet mode (no output)\n");
	exit(0);
}

int main(int argc, char **argv)
{
	long *predy;
	double performance;

	char input_filename[MAX_LINE_LENGTH];
	char model_filename[MAX_LINE_LENGTH];
	char output_filename[MAX_LINE_LENGTH];;

	struct MajModel *model = Malloc(struct MajModel, 1);
	struct MajData *data = Malloc(struct MajData, 1);

	if (argc < MINARGS || msvmmaj_check_argv(argc, argv, "-help") 
			|| msvmmaj_check_argv_eq(argc, argv, "-h") )
		exit_with_help();
	parse_command_line(argc, argv, model, input_filename, output_filename,
			model_filename);

	// TODO: make sure that read_data allows for files without labels
	msvmmaj_read_data(data, input_filename);
	msvmmaj_read_model(model, model_filename);

	// check if the number of attributes in data equals that in model
	if (data->m != model->m) {
		fprintf(stderr, "Error: number of attributes in data (%li) "
				"does not equal the number of attributes in "
				"model (%li)\n", data->m, model->m);
		exit(1);
	}

	predy = Calloc(long, data->n);
	msvmmaj_predict_labels(data, model, predy);
	if (data->y != NULL) {
		performance = msvmmaj_prediction_perf(data, predy);
		note("Predictive performance: %3.2f%%\n", performance);
	}

	if (msvmmaj_check_argv_eq(argc, argv, "-o")) {
		msvmmaj_write_predictions(data, predy, output_filename);
		note("Predictions written to: %s\n", output_filename);
	}

	msvmmaj_free_model(model);
	msvmmaj_free_data(data);
	free(predy);

	return 0;
}

void parse_command_line(int argc, char **argv, struct MajModel *model, 
		char *input_filename, char *output_filename, char *model_filename)
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
				fprintf(stderr, "Unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
		}
	}

	if (i >= argc)
		exit_with_help();

	strcpy(input_filename, argv[i]);
	i++;
	strcpy(model_filename, argv[i]);
}
