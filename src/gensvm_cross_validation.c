/**
 * @file gensvm_cross_validation.c
 * @author G.J.J. van den Burg
 * @date 2016-10-24
 * @brief Function for running cross validation on GenModel
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

#include "gensvm_cross_validation.h"

extern FILE *GENSVM_OUTPUT_FILE;

/**
 * @brief Run cross validation with a given set of train/test folds
 *
 * @details
 * This cross validation function uses predefined train/test splits. Also, the
 * the optimal parameters GenModel::V of a previous fold as initial conditions
 * for GenModel::V of the next fold.
 *
 * @note
 * This function always sets the output stream defined in GENSVM_OUTPUT_FILE
 * to NULL, to ensure gensvm_optimize() doesn't print too much.
 *
 * @param[in] 	model 		GenModel with the configuration to train
 * @param[in] 	train_folds 	array of training datasets
 * @param[in] 	test_folds 	array of test datasets
 * @param[in] 	folds 		number of folds
 * @param[in] 	n_total 	number of objects in the union of the train
 * 				datasets
 * @return 			performance (hitrate) of the configuration on
 * 				cross validation
 */
double gensvm_cross_validation(struct GenModel *model,
		struct GenData **train_folds, struct GenData **test_folds,
		long folds, long n_total)
{
	long f;
	long *predy = NULL;
	double performance, total_perf = 0;

	// make sure that gensvm_optimize() is silent.
	FILE *fid = GENSVM_OUTPUT_FILE;
	GENSVM_OUTPUT_FILE = NULL;

	// run cross-validation
	for (f=0; f<folds; f++) {
		// reallocate model in case dimensions differ with data
		gensvm_reallocate_model(model, train_folds[f]->n,
				train_folds[f]->r);

		// initialize object weights
		gensvm_initialize_weights(train_folds[f], model);

		// train the model (surpressing output)
		gensvm_optimize(model, train_folds[f]);

		// calculate prediction performance on test set
		predy = Calloc(long, test_folds[f]->n);
		gensvm_predict_labels(test_folds[f], model, predy);
		performance = gensvm_prediction_perf(test_folds[f], predy);
		total_perf += performance * test_folds[f]->n;

		free(predy);
	}

	total_perf /= ((double) n_total);

	// reset the output stream
	GENSVM_OUTPUT_FILE = fid;

	return total_perf;
}

void copy_predictions(long *predy, long *predictions, long *cv_idx, long fold, 
		long N)
{
	long i, j = 0;
	for (i=0; i<N; i++) {
		if (cv_idx[i] == fold) {
			predictions[i] = predy[j];
			j++;
		}
	}
}

void gensvm_cross_validation_store(struct GenModel *model, 
		struct GenData **train_folds, struct GenData **test_folds, 
		long folds, long n_total, long *cv_idx, long *predictions, 
		double *durations, int verbosity)
{
	FILE *fid = NULL;
	long f;
	long *predy = NULL;
	void (*tmpfunc) (const char *, ...) = NULL;

	GenTime fold_s, fold_e;

	gensvm_py_reset_interrupt_hdl();

	// make sure that gensvm_optimize() is silent.
	if (verbosity <= 1) {
		fid = GENSVM_OUTPUT_FILE;
		GENSVM_OUTPUT_FILE = NULL;
		tmpfunc = gensvm_print_out;
		gensvm_print_out = gensvm_print_output_fpt;
	}

	// run cross-validation
	for (f=0; f<folds; f++) {
		Timer(fold_s);

		// reallocate model in case dimensions differ with data
		gensvm_reallocate_model(model, train_folds[f]->n,
				train_folds[f]->r);

		// initialize object weights
		gensvm_initialize_weights(train_folds[f], model);

		// train the model (surpressing output)
		gensvm_optimize(model, train_folds[f]);

		// calculate prediction performance on test set
		predy = Calloc(long, test_folds[f]->n);
		gensvm_predict_labels(test_folds[f], model, predy);
		copy_predictions(predy, predictions, cv_idx, f, n_total);
		free(predy);

		Timer(fold_e);
		durations[f] = gensvm_elapsed_time(&fold_s, &fold_e);

		if (gensvm_py_pending_interrupt()) {
			break;
		}
	}

	// reset the output stream
	if (verbosity <= 1) {
		GENSVM_OUTPUT_FILE = fid;
		gensvm_print_out = tmpfunc;
	}
}
