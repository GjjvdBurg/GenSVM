/**
 * @file gensvm_gridsearch.c
 * @author G.J.J. van den Burg
 * @date 2014-01-07
 * @brief Functions for finding the optimal parameters for the dataset
 *
 * @details
 * The GenSVM algorithm takes a number of parameters. The functions in
 * this file are used to find the optimal parameters.
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

#include "gensvm_gridsearch.h"

extern FILE *GENSVM_OUTPUT_FILE;

/**
 * @brief Initialize a GenQueue from a Training instance
 *
 * @details
 * A Training instance describes the grid to search over. This funtion
 * creates all tasks that need to be performed and adds these to
 * a GenQueue. Each task contains a pointer to the train and test datasets
 * which are supplied. Note that the tasks are created in a specific order of
 * the parameters, to ensure that the GenModel::V of a previous parameter
 * set provides the best possible initial estimate of GenModel::V for the next
 * parameter set.
 *
 * @param[in] 	grid 	Training struct describing the grid search
 * @param[in] 	queue 		pointer to a GenQueue that will be used to
 * 				add the tasks to
 * @param[in] 	train_data 	GenData of the training set
 * @param[in] 	test_data 	GenData of the test set
 *
 */
void gensvm_fill_queue(struct GenGrid *grid, struct GenQueue *queue,
		struct GenData *train_data, struct GenData *test_data)
{
	long i, j, k;
	long N, cnt = 0;
	struct GenTask *task = NULL;
	queue->i = 0;

	N = grid->Np;
	N *= grid->Nl;
	N *= grid->Nk;
	N *= grid->Ne;
	N *= grid->Nw;
	// these parameters are not necessarily non-zero
	N *= grid->Ng > 0 ? grid->Ng : 1;
	N *= grid->Nc > 0 ? grid->Nc : 1;
	N *= grid->Nd > 0 ? grid->Nd : 1;

	queue->tasks = Calloc(struct GenTask *, N);
	queue->N = N;

	// initialize all tasks
	for (i=0; i<N; i++) {
		task = gensvm_init_task();
		task->ID = i;
		task->train_data = train_data;
		task->test_data = test_data;
		task->folds = grid->folds;
		task->kerneltype = grid->kerneltype;
		queue->tasks[i] = task;
	}

	// These loops mimick a large nested for loop. The advantage is that
	// Nd, Nc and Ng which are on the outside of the nested for loop can
	// now be zero, without large modification (see below). Whether this
	// is indeed better than the nested for loop has not been tested.
	cnt = 1;
	i = 0;
	while (i < N) {
		for (j=0; j<grid->Np; j++) {
			for (k=0; k<cnt; k++) {
				queue->tasks[i]->p = grid->ps[j];
				i++;
			}
		}
	}

	cnt *= grid->Np;
	i = 0;
	while (i < N) {
		for (j=0; j<grid->Nl; j++) {
			for (k=0; k<cnt; k++) {
				queue->tasks[i]->lambda =
					grid->lambdas[j];
				i++;
			}
		}
	}

	cnt *= grid->Nl;
	i = 0;
	while (i < N) {
		for (j=0; j<grid->Nk; j++) {
			for (k=0; k<cnt; k++) {
				queue->tasks[i]->kappa = grid->kappas[j];
				i++;
			}
		}
	}

	cnt *= grid->Nk;
	i = 0;
	while (i < N) {
		for (j=0; j<grid->Nw; j++) {
			for (k=0; k<cnt; k++) {
				queue->tasks[i]->weight_idx =
					grid->weight_idxs[j];
				i++;
			}
		}
	}

	cnt *= grid->Nw;
	i = 0;
	while (i < N) {
		for (j=0; j<grid->Ne; j++) {
			for (k=0; k<cnt; k++) {
				queue->tasks[i]->epsilon = grid->epsilons[j];
				i++;
			}
		}
	}

	cnt *= grid->Ne;
	i = 0;
	while (i < N && grid->Ng > 0) {
		for (j=0; j<grid->Ng; j++) {
			for (k=0; k<cnt; k++) {
				queue->tasks[i]->gamma = grid->gammas[j];
				i++;
			}
		}
	}

	cnt *= grid->Ng > 0 ? grid->Ng : 1;
	i = 0;
	while (i < N && grid->Nc > 0) {
		for (j=0; j<grid->Nc; j++) {
			for (k=0; k<cnt; k++) {
				queue->tasks[i]->coef = grid->coefs[j];
				i++;
			}
		}
	}

	cnt *= grid->Nc > 0 ? grid->Nc : 1;
	i = 0;
	while (i < N && grid->Nd > 0) {
		for (j=0; j<grid->Nd; j++) {
			for (k=0; k<cnt; k++) {
				queue->tasks[i]->degree = grid->degrees[j];
				i++;
			}
		}
	}
}

/**
 * @brief Check if the kernel parameters change between tasks
 *
 * @details
 * In the current strategy for training the kernel matrix is decomposed once,
 * and tasks with the same kernel settings are performed sequentially. When a
 * task needs to be done with different kernel parameters, the kernel matrix
 * needs to be recalculated. This function is used to check whether this is
 * the case.
 *
 * @param[in] 	newtask 	the next task
 * @param[in] 	oldtask 	the old task
 * @return 			whether the kernel needs to be reevaluated
 */
bool gensvm_kernel_changed(struct GenTask *newtask, struct GenTask *oldtask)
{
	if (oldtask == NULL)
		return true;
	if (newtask->kerneltype != oldtask->kerneltype) {
		return true;
	} else if (newtask->kerneltype == K_POLY) {
		if (newtask->gamma != oldtask->gamma)
			return true;
		if (newtask->coef != oldtask->coef)
			return true;
		if (newtask->degree != oldtask->degree)
			return true;
		return false;
	} else if (newtask->kerneltype == K_RBF) {
		if (newtask->gamma != oldtask->gamma)
			return true;
		return false;
	} else if (newtask->kerneltype == K_SIGMOID) {
		if (newtask->gamma != oldtask->gamma)
			return true;
		if (newtask->coef != oldtask->coef)
			return true;
		return false;
	}
	return false;
}

/**
 * @brief Compute the kernels for the folds of the train and test datasets
 *
 * @details
 * When the kernel parameters change in a kernel grid search, the kernel
 * pre- and post-processing has to be done for the new kernel parameters. This
 * is done here for each of the folds. Each of the training folds is
 * preprocessed, and each of the test folds is postprocessed.
 *
 * @param[in] 		folds 		number of cross validation folds
 * @param[in] 		model 		GenModel with new kernel parameters
 * @param[in,out] 	train_folds 	array of train datasets
 * @param[in,out] 	test_folds 	array of test datasets
 *
 */
void gensvm_kernel_folds(long folds, struct GenModel *model,
		struct GenData **train_folds, struct GenData **test_folds)
{
	long f;

	if (model->kerneltype != K_LINEAR)
		note("Computing kernel ... ");
	for (f=0; f<folds; f++) {
		if (train_folds[f]->Z != train_folds[f]->RAW)
			free(train_folds[f]->Z);
		if (test_folds[f]->Z != test_folds[f]->RAW)
			free(test_folds[f]->Z);
		gensvm_kernel_preprocess(model, train_folds[f]);
		gensvm_kernel_postprocess(model, train_folds[f],
				test_folds[f]);
	}
	if (model->kerneltype != K_LINEAR)
		note("done.\n");
}

/**
 * @brief Run the grid search for a GenQueue
 *
 * @details
 * Given a GenQueue of GenTask struct to be trained, a grid search is launched to
 * find the optimal parameter configuration. As is also done within
 * cross_validation(), the optimal weights of one parameter set are used as
 * initial estimates for GenModel::V in the next parameter set. Note that to
 * optimally exploit this feature of the optimization algorithm, the order in
 * which tasks are considered is important. This is considered in
 * make_queue().
 *
 * The performance found by cross validation is stored in the GenTask struct.
 *
 * @param[in,out] 	q 			GenQueue with GenTask
 * 						instances to run
 * @param[in] 		cv_idx 			cross validation index array
 * 						(optional)
 * @param[in] 		store_predictions 	whether or not to store the
 * 						full array of cross validation
 * 						predictions
 * @param[in] 		verbose 		set the verbosity level.
 *
 * @return 	total training time
 */
double gensvm_train_queue(struct GenQueue *q, long *cv_idx,
		bool store_predictions, int verbose)
{
	bool free_cv_idx = false;
	long f, folds;
	double perf, total_time, current_max = -1;
	struct GenTask *task = get_next_task(q);
	struct GenTask *prevtask = NULL;
	struct GenModel *model = gensvm_init_model();
	GenTime main_s, main_e, loop_s, loop_e;

	folds = task->folds;

	gensvm_py_reset_interrupt_hdl();

	model->n = 0;
	model->m = task->train_data->m;
	model->K = task->train_data->K;
	gensvm_allocate_model(model);
	gensvm_init_V(NULL, model, task->train_data);

	if (cv_idx == NULL) {
		cv_idx = Calloc(long, task->train_data->n);
		gensvm_make_cv_split(task->train_data->n, folds, cv_idx);
		free_cv_idx = true;
	}

	struct GenData **train_folds = Malloc(struct GenData *, folds);
	struct GenData **test_folds = Malloc(struct GenData *, folds);
	for (f=0; f<folds; f++) {
		train_folds[f] = gensvm_init_data();
		test_folds[f] = gensvm_init_data();
		gensvm_get_tt_split(task->train_data, train_folds[f],
				test_folds[f], cv_idx, f);
	}

	note("Starting grid search ...\n");
	Timer(main_s);
	while (task) {
		gensvm_task_to_model(task, model);
		if (gensvm_kernel_changed(task, prevtask)) {
			gensvm_kernel_folds(folds, model, train_folds,
					test_folds);
		}

		if (store_predictions) {
			long *predictions = Calloc(long, task->train_data->n);
			for (f=0; f<task->train_data->n; f++)
				predictions[f] = -1;
			double *durations = Calloc(double, folds);
			for (f=0; f<folds; f++)
				durations[f] = -1;
			Timer(loop_s);
			gensvm_cross_validation_store(model, train_folds,
					test_folds, folds, task->train_data->n,
					cv_idx, predictions, durations,
					verbose);
			Timer(loop_e);
			task->predictions = predictions;
			task->durations = durations;
		}
		else {
			Timer(loop_s);
			perf = gensvm_cross_validation(model, train_folds,
					test_folds, folds,
					task->train_data->n);
			Timer(loop_e);
			current_max = maximum(current_max, perf);
			task->performance = perf;
		}

		task->duration = gensvm_elapsed_time(&loop_s, &loop_e);

		gensvm_gridsearch_progress(task, q->N, task->performance,
				task->duration, current_max, !store_predictions);

		prevtask = task;
		task = get_next_task(q);

		if (gensvm_py_pending_interrupt())
			break;
	}
	Timer(main_e);

	total_time = gensvm_elapsed_time(&main_s, &main_e);
	note("\nTotal time: %8.8f seconds\n", total_time);

	gensvm_free_model(model);
	for (f=0; f<folds; f++) {
		gensvm_free_data(train_folds[f]);
		gensvm_free_data(test_folds[f]);
	}
	free(train_folds);
	free(test_folds);
	if (free_cv_idx)
		free(cv_idx);

	return total_time;
}

/**
 * @brief Print the description of the current task on screen
 *
 * @details
 * To track the progress of the grid search the parameters of the current task
 * are written to the output specified in GENSVM_OUTPUT_FILE. Since the
 * parameters differ with the specified kernel, this function writes a
 * parameter string depending on which kernel is used.
 *
 * @param[in] 	task 		the GenTask specified
 * @param[in] 	N 		total number of tasks
 * @param[in] 	perf 		performance of the current task
 * @param[in] 	duration 	time duration of the current task
 * @param[in] 	current_max 	current best performance
 *
 */
void gensvm_gridsearch_progress(struct GenTask *task, long N, double perf,
		double duration, double current_max, bool show_perf)
{
	char buffer[GENSVM_MAX_LINE_LENGTH];
	sprintf(buffer, "(%03li/%03li)\t", task->ID+1, N);
	if (task->kerneltype == K_POLY)
		sprintf(buffer + strlen(buffer), "d = %2.2f\t", task->degree);
	if (task->kerneltype == K_POLY || task->kerneltype == K_SIGMOID)
		sprintf(buffer + strlen(buffer), "c = %2.2f\t", task->coef);
	if (task->kerneltype == K_POLY || task->kerneltype == K_SIGMOID ||
			task->kerneltype == K_RBF)
		sprintf(buffer + strlen(buffer), "g = %3.3f\t", task->gamma);
	sprintf(buffer + strlen(buffer), "eps = %g\tw = %i\tk = %2.2f\t"
			"l = %11g\tp = %2.2f\t", task->epsilon,
			task->weight_idx, task->kappa, task->lambda, task->p);
	note(buffer);
	if (show_perf)
		note("%3.3f%% (%3.3fs)\t(best = %3.3f%%)\n", perf, duration,
			current_max);
	else
		note("(%3.3fs)\n", duration);
}
