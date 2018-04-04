/**
 * @file gensvm_consistency.c
 * @author G.J.J. van den Burg
 * @date 2016-10-24
 * @brief Functions for running consistency repeats
 *
 * @details
 * When running the grid search, the user may optionally choose to run 
 * repetitions of the best performing configurations to check if they perform 
 * well consistently. These are called consistency repeats, and the code 
 * needed for them is defined here.
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

#include "gensvm_consistency.h"

/**
 * @brief Create GenQueue of tasks with performance above a given percentile
 *
 * @details
 * This function constructs a GenQueue of the GenTask instances in the 
 * provided GenQueue which have a performance at or above the given percentile 
 * of the performance of all elements in the provided GenQueue. This can be 
 * used to determine which hyperparameter configurations belong to the top x-% 
 * of all tasks in terms of performance.
 *
 * @sa
 * gensvm_consistency_repeats(), gensvm_percentile()
 *
 * @note
 * This function assumes that for each task in the given GenQueue, the 
 * GenTask::perf element has been set.
 *
 * @param[in] 	q 		a complete GenQueue struct
 * @param[in] 	percentile 	the desired percentile
 *
 * @return 			a GenQueue struct with GenTasks which are at
 * 				or above the desired percentile of performance
 *
 */
struct GenQueue *gensvm_top_queue(struct GenQueue *q, double percentile)
{
	long i, k, N = 0;
	double boundary,
	       *perf = Calloc(double, q->N);
	struct GenQueue *nq = gensvm_init_queue();

	// find the desired percentile of performance
	for (i=0; i<q->N; i++) {
		perf[i] = q->tasks[i]->performance;
	}
	boundary = gensvm_percentile(perf, q->N, percentile);
	note("Boundary of the %g-th percentile determined at: %f\n",
			percentile, boundary);

	// find the number of tasks that perform at or above the boundary
	for (i=0; i<q->N; i++) {
		if (q->tasks[i]->performance >= boundary)
			N++;
	}

	// create a new queue with the best tasks
	nq->tasks = Malloc(struct GenTask *, N);
	k = 0;
	for (i=0; i<q->N; i++) {
		if (q->tasks[i]->performance >= boundary)
			nq->tasks[k++] = gensvm_copy_task(q->tasks[i]);
	}
	nq->N = N;
	nq->i = 0;

	free(perf);
	return nq;
}

/**
 * @brief Run repeats of the GenTask structs in GenQueue to find the best
 * configuration
 *
 * @details
 * The best performing tasks in the supplied GenQueue are found by taking 
 * those GenTask structs that have a performance greater or equal to the given 
 * percentile of the performance of all tasks. These tasks are then gathered 
 * in a new GenQueue. For each of the tasks in this new GenQueue the cross 
 * validation run is repeated a number of times.
 *
 * For each of the GenTask configurations that are repeated the mean 
 * performance, standard deviation of the performance and the mean computation 
 * time are reported.
 *
 * Finally, the overall best tasks are written to the specified output. These 
 * tasks are selected to have both the highest mean performance, as well as 
 * the smallest standard deviation in their performance. This is done as 
 * follows.  First the 99th percentile of task performance and the 1st 
 * percentile of standard deviation is calculated. If a task exists for which 
 * the mean performance of the repeats and the standard deviation equals these 
 * values respectively, this task is found to be the best and is written to 
 * the output. If no such task exists, the 98th percentile of performance and 
 * the 2nd percentile of standard deviation is considered. This is repeated 
 * until an interval is found which contains tasks. If one or more tasks are 
 * found, this loop stops.
 *
 * @param[in] 	q 		GenQueue of GenTask structs which have already been
 * 				run and have a GenTask::performance value
 * @param[in] 	repeats 	Number of times to repeat the best
 * 				configurations for consistency
 * @param[in] 	percentile 	percentile of performance to determine which
 * 				tasks to repeat
 *
 * @return 			ID of the best task
 *
 */
int gensvm_consistency_repeats(struct GenQueue *q, long repeats,
		double percentile)
{
	bool breakout;
	long i, f, r, N, *cv_idx = NULL;
	double p, pi, pr, pt,
	       *time = NULL,
	       *std = NULL,
	       *mean = NULL,
	       *perf = NULL;
	struct GenQueue *nq = NULL;
	struct GenData **train_folds = NULL,
		       **test_folds = NULL;
	struct GenModel *model = gensvm_init_model();
	struct GenTask *task = NULL;
	GenTime loop_s, loop_e;

	nq = gensvm_top_queue(q, percentile);
	N = nq->N;

	note("Number of items to check: %li\n", nq->N);
	std = Calloc(double, N);
	mean = Calloc(double, N);
	time = Calloc(double, N);
	perf = Calloc(double, N*repeats);

	task = get_next_task(nq);

	model->n = 0;
	model->m = task->train_data->m;
	model->K = task->train_data->K;
	gensvm_allocate_model(model);
	gensvm_init_V(NULL, model, task->train_data);

	cv_idx = Calloc(long, task->train_data->n);

	train_folds = Malloc(struct GenData *, task->folds);
	test_folds = Malloc(struct GenData *, task->folds);

	i = 0;
	while (task) {
		gensvm_task_to_model(task, model);

		time[i] = 0.0;
		note("(%02li/%02li:%03li)\t", i+1, N, task->ID);
		for (r=0; r<repeats; r++) {
			Memset(cv_idx, long, task->train_data->n);
			gensvm_make_cv_split(task->train_data->n, task->folds, cv_idx);
			train_folds = Malloc(struct GenData *, task->folds);
			test_folds = Malloc(struct GenData *, task->folds);
			for (f=0; f<task->folds; f++) {
				train_folds[f] = gensvm_init_data();
				test_folds[f] = gensvm_init_data();
				gensvm_get_tt_split(task->train_data, train_folds[f],
						test_folds[f], cv_idx, f);
				gensvm_kernel_preprocess(model, train_folds[f]);
				gensvm_kernel_postprocess(model, train_folds[f],
						test_folds[f]);
			}

			Timer(loop_s);
			p = gensvm_cross_validation(model, train_folds, test_folds,
					task->folds, task->train_data->n);
			Timer(loop_e);
			time[i] += gensvm_elapsed_time(&loop_s, &loop_e);
			matrix_set(perf, N, repeats, i, r, p);
			mean[i] += p/((double) repeats);
			note("%3.3f\t", p);
			// this is done because if we reuse the V it's not a
			// consistency check
			gensvm_init_V(NULL, model, task->train_data);
			for (f=0; f<task->folds; f++) {
				gensvm_free_data(train_folds[f]);
				gensvm_free_data(test_folds[f]);
			}
			free(train_folds);
			train_folds = NULL;

			free(test_folds);
			test_folds = NULL;
		}
		for (r=0; r<repeats; r++) {
			std[i] += pow(matrix_get(perf, N, repeats, i, r) 
					- mean[i], 2.0);
		}
		if (r > 1) {
			std[i] /= ((double) repeats) - 1.0;
			std[i] = sqrt(std[i]);
		} else {
			std[i] = 0.0;
		}
		note("(m = %3.3f, s = %3.3f, t = %3.3f)\n", mean[i], std[i],
				time[i]);
		task = get_next_task(nq);
		i++;
	}

	// find the best overall configurations: those with high average
	// performance and low deviation in the performance
	note("\nBest overall configuration(s):\n");
	note("ID\tweights\tepsilon\t\tp\t\tkappa\t\tlambda\t\t"
			"mean_perf\tstd_perf\ttime_perf\n");
	p = 0.0;
	breakout = false;
	int best_id = -1;
	while (!breakout) {
		pi = gensvm_percentile(mean, N, (100.0-p));
		pr = gensvm_percentile(std, N, p);
		pt = gensvm_percentile(time, N, p);
		for (i=0; i<N; i++) {
			if ((pi - mean[i] < 0.0001) &&
					(std[i] - pr < 0.0001) &&
					(time[i] - pt < 0.0001)) {
				note("(%li)\tw = %li\te = %f\tp = %f\t"
						"k = %f\tl = %f\t"
						"mean: %3.3f\tstd: %3.3f\t"
						"time: %3.3f\n",
						nq->tasks[i]->ID,
						nq->tasks[i]->weight_idx,
						nq->tasks[i]->epsilon,
						nq->tasks[i]->p,
						nq->tasks[i]->kappa,
						nq->tasks[i]->lambda,
						mean[i],
						std[i],
						time[i]);
				breakout = true;
				if (best_id == -1)
					best_id = nq->tasks[i]->ID;
			}
		}
		p += 1.0;
	}

	free(cv_idx);
	gensvm_free_model(model);
	gensvm_free_queue(nq);

	free(perf);
	free(std);
	free(mean);
	free(time);

	return best_id;
}

/**
 * @brief Comparison function for doubl
 *
 * @param[in] 	elem1 	number 1
 * @param[in] 	elem2 	number 2
 * @returns 		comparison of number 1 larger than number 2
 */
int gensvm_dsort(const void *elem1, const void *elem2)
{
	const double t1 = (*(double *) elem1);
	const double t2 = (*(double *) elem2);
	return t1 > t2;
}

/**
 * @brief Calculate the percentile of an array of doubles
 *
 * @details
 * The percentile of performance is used to find the top performing
 * configurations. Since no standard definition of the percentile exists, we
 * use the method used in MATLAB and Octave. Since calculating the percentile
 * requires a sorted list of the values, a local copy is made first.
 *
 * @param[in] 	values 	array of doubles
 * @param[in] 	N 	length of the array
 * @param[in] 	p 	percentile to calculate ( 0 <= p <= 100.0 ).
 * @returns 		the p-th percentile of the values
 */
double gensvm_percentile(double *values, long N, double p)
{
	if (N == 1)
		return values[0];

	long i;
	double pi, pr, boundary;
	double *local = Malloc(double, N);
	for (i=0; i<N; i++)
		local[i] = values[i];

	qsort(local, N, sizeof(double), gensvm_dsort);
	p /= 100.0;
	p = p*N + 0.5;
	pi = maximum(minimum(floor(p), N-1), 1);
	pr = maximum(minimum(p - pi, 1), 0);
	boundary = (1 - pr)*local[((long) pi)-1] + pr*local[((long) pi)];

	free(local);

	return boundary;
}
