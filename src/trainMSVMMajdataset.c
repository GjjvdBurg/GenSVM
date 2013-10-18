#include <time.h>

#include "crossval.h"
#include "MSVMMaj.h"
#include "msvmmaj_pred.h"
#include "msvmmaj_train.h"
#include "msvmmaj_train_dataset.h"
#include "strutil.h"
#include "util.h"

#define MINARGS 2

extern FILE *MSVMMAJ_OUTPUT_FILE;

void print_null(const char *s) {}
void exit_with_help();
void parse_command_line(int argc, char **argv, char *input_filename);
void read_training_from_file(char *input_filename, struct Training *training);

void exit_with_help()
{
	printf("This is MSVMMaj, version %1.1f\n\n", VERSION);
	printf("Usage: trainMSVMMajdataset [options] training_file\n");
	printf("Options:\n");
	printf("-h | -help : print this help.\n");
	printf("-q : quiet mode (no output)\n");

	exit(0);
}

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

	free_queue(q);
	free(training);
	msvmmaj_free_data(train_data);
	msvmmaj_free_data(test_data);

	note("Done.\n");
	return 0;
}

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
				fprintf(stderr, "Unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
		}
	}

	if (i >= argc)
		exit_with_help();

	strcpy(input_filename, argv[i]);
}

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
		fprintf(stderr, "Error opening training file %s\n", input_filename);
		exit(1);
	}
	training->traintype = CV;
	while ( fgets(buffer, MAX_LINE_LENGTH, fid) != NULL ) {
		Memset(params, double,  MAX_LINE_LENGTH);
		Memset(lparams, long, MAX_LINE_LENGTH);
		if (str_startswith(buffer, "train:")) {
			sscanf(buffer, "train: %s\n", train_filename);
			training->train_data_file = Calloc(char, MAX_LINE_LENGTH);
			strcpy(training->train_data_file, train_filename);
		} else if (str_startswith(buffer, "test:")) {
			sscanf(buffer, "test: %s\n", test_filename);
			training->test_data_file = Calloc(char, MAX_LINE_LENGTH);
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
				fprintf(stderr, "Field \"folds\" only takes one value. "
						"Additional fields are ignored.\n");
		} else if (str_startswith(buffer, "repeats:")) {
			nr = all_longs_str(buffer, 8, lparams);
			training->repeats = lparams[0];
			if (nr > 1)
				fprintf(stderr, "Field \"repeats\" only takes one value. "
						"Additional fields are ignored.\n");
		} else {
			fprintf(stderr, "Cannot find any parameters on line: %s\n", buffer);
		}
	}

	free(params);
	free(lparams);
	fclose(fid);
}
