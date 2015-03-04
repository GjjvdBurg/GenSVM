/**
 * @file gensvm_train_dataset.c
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Functions for finding the optimal parameters for the dataset
 *
 * @details
 * The GenSVM algorithm takes a number of parameters. The functions in
 * this file are used to find the optimal parameters. 
 */

#include <math.h>
#include <time.h>

#include "crossval.h"
#include "libGenSVM.h"
#include "gensvm.h"
#include "gensvm_init.h"
#include "gensvm_kernel.h"
#include "gensvm_matrix.h"
#include "gensvm_train.h"
#include "gensvm_train_dataset.h"
#include "gensvm_pred.h"
#include "util.h"
#include "timer.h"

extern FILE *GENSVM_OUTPUT_FILE;

/**
 * @brief Initialize a Queue from a Training instance
 *
 * @details
 * A Training instance describes the grid to search over. This funtion
 * creates all tasks that need to be performed and adds these to 
 * a Queue. Each task contains a pointer to the train and test datasets
 * which are supplied. Note that the tasks are created in a specific order of
 * the parameters, to ensure that the GenModel::V of a previous parameter
 * set provides the best possible initial estimate of GenModel::V for the next
 * parameter set.
 *
 * @param[in] 	training 	Training struct describing the grid search
 * @param[in] 	queue 		pointer to a Queue that will be used to 
 * 				add the tasks to
 * @param[in] 	train_data 	GenData of the training set
 * @param[in] 	test_data 	GenData of the test set
 *
 */
void make_queue(struct Training *training, struct Queue *queue, 
		struct GenData *train_data, struct GenData *test_data)
{
	long i, j, k;
	long N, cnt = 0;
	struct Task *task;
	queue->i = 0;

	N = training->Np;
	N *= training->Nl;
	N *= training->Nk;
	N *= training->Ne;
	N *= training->Nw;
	// these parameters are not necessarily non-zero
	N *= training->Ng > 0 ? training->Ng : 1;
	N *= training->Nc > 0 ? training->Nc : 1;
	N *= training->Nd > 0 ? training->Nd : 1;

	queue->tasks = Calloc(struct Task *, N);
	queue->N = N;

	// initialize all tasks
	for (i=0; i<N; i++) {
		task = Calloc(struct Task, 1);
		task->ID = i;
		task->train_data = train_data;
		task->test_data = test_data;
		task->folds = training->folds;
		task->kerneltype = training->kerneltype;
		task->kernelparam = Calloc(double, training->Ng + 
				training->Nc + training->Nd);
		queue->tasks[i] = task;
	}

	// These loops mimick a large nested for loop. The advantage is that
	// Nd, Nc and Ng which are on the outside of the nested for loop can
	// now be zero, without large modification (see below). Whether this
	// is indeed better than the nested for loop has not been tested.
	cnt = 1;
	i = 0;
	while (i < N )
		for (j=0; j<training->Np; j++) 
			for (k=0; k<cnt; k++) {
				queue->tasks[i]->p = training->ps[j];
				i++;
			}

	cnt *= training->Np;
	i = 0;
	while (i < N )
		for (j=0; j<training->Nl; j++) 
			for (k=0; k<cnt; k++) {
				queue->tasks[i]->lambda = 
					training->lambdas[j];
				i++;
			}

	cnt *= training->Nl;
	i = 0;
	while (i < N )
		for (j=0; j<training->Nk; j++)
			for (k=0; k<cnt; k++) {
				queue->tasks[i]->kappa = training->kappas[j];
				i++;
			}

	cnt *= training->Nk;
	i = 0;
	while (i < N )
		for (j=0; j<training->Nw; j++)
			for (k=0; k<cnt; k++) {
				queue->tasks[i]->weight_idx = 
					training->weight_idxs[j];
				i++;
			}

	cnt *= training->Nw;
	i = 0;
	while (i < N )
		for (j=0; j<training->Ne; j++)
			for (k=0; k<cnt; k++) {
				queue->tasks[i]->epsilon = 
					training->epsilons[j];
				i++;
			}

	cnt *= training->Ne;
	i = 0;
	while (i < N && training->Ng > 0)
		for (j=0; j<training->Ng; j++)
			for (k=0; k<cnt; k++) {
				queue->tasks[i]->kernelparam[0] = 
					training->gammas[j];
				i++;
			}

	cnt *= training->Ng > 0 ? training->Ng : 1;
	i = 0;
	while (i < N && training->Nc > 0)
		for (j=0; j<training->Nc; j++)
			for (k=0; k<cnt; k++) {
				queue->tasks[i]->kernelparam[1] = 
					training->coefs[j];
				i++;
			}

	cnt *= training->Nc > 0 ? training->Nc : 1;
	i = 0;
	while (i < N && training->Nd > 0)
		for (j=0; j<training->Nd; j++)
			for (k=0; k<cnt; k++) {
				queue->tasks[i]->kernelparam[2] = 
					training->degrees[j];
				i++;
			}
}

/**
 * @brief Get new Task from Queue
 *
 * @details
 * Return a pointer to the next Task in the Queue. If no Task instances are
 * left, NULL is returned. The internal counter Queue::i is used for finding
 * the next Task.
 *
 * @param[in] 	q 	Queue instance
 * @returns 		pointer to next Task
 *
 */
struct Task *get_next_task(struct Queue *q)
{
	long i = q->i;
	if (i < q->N) {
		q->i++;
		return q->tasks[i];
	}
	return NULL;
}

/**
 * @brief Comparison function for Tasks based on performance
 *
 * @details
 * To be able to sort Task structures on the performance of their specific
 * set of parameters, this comparison function is implemented. Task structs
 * are sorted with highest performance first.
 *
 * @param[in] 	elem1 	Task 1
 * @param[in] 	elem2 	Task 2
 * @returns 		result of inequality of Task 1 performance over
 * 			Task 2 performance
 */
int tasksort(const void *elem1, const void *elem2)
{
	const struct Task *t1 = (*(struct Task **) elem1);
	const struct Task *t2 = (*(struct Task **) elem2);
	return (t1->performance > t2->performance);
}

/**
 * @brief Comparison function for doubl
 *
 * @param[in] 	elem1 	number 1
 * @param[in] 	elem2 	number 2
 * @returns 		comparison of number 1 larger than number 2
 */
int doublesort(const void *elem1, const void *elem2)
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
double prctile(double *values, long N, double p)
{
	if (N == 1)
		return values[0];

	long i;
	double pi, pr, boundary;
	double *local = Malloc(double, N);
	for (i=0; i<N; i++)
		local[i] = values[i];

	qsort(local, N, sizeof(double), doublesort);
	p /= 100.0;
	p = p*N + 0.5;
	pi = maximum(minimum(floor(p), N-1), 1);
	pr = maximum(minimum(p - pi, 1), 0);
	boundary = (1 - pr)*local[((long) pi)-1] + pr*local[((long) pi)];

	free(local);

	return boundary;
}

struct Queue *create_top_queue(struct Queue *q)
{
	long i, k, N = 0;
	double boundary, *perf;
	struct Queue *nq = Malloc(struct Queue, 1);

	// find the 95th percentile of performance
	perf = Calloc(double, q->N);
	for (i=0; i<q->N; i++) {
		perf[i] = q->tasks[i]->performance;
	}
	boundary = prctile(perf, q->N, 95.0);
	free(perf);
	note("boundary determined at: %f\n", boundary);

	// find the number of tasks that perform at or above the boundary
	for (i=0; i<q->N; i++) {
		if (q->tasks[i]->performance >= boundary)
			N++;
	}

	// create a new queue with the best tasks
	nq->tasks = Malloc(struct Task *, N);
	k = 0;
	for (i=0; i<q->N; i++) {
		if (q->tasks[i]->performance >= boundary)
			nq->tasks[k++] = q->tasks[i];
	}
	nq->N = N;
	nq->i = 0;
	
	return nq;
}


/**
 * @brief Run repeats of the Task structs in Queue to find the best
 * configuration
 *
 * @details
 * The best performing tasks in the supplied Queue are found by taking those
 * Task structs that have a performance greater or equal to the 95% percentile
 * of the performance of all tasks. These tasks are then gathered in a new
 * Queue. For each of the tasks in this new Queue the cross validation run is 
 * repeated a number of times. 
 *
 * For each of the Task configurations that are repeated the mean performance,
 * standard deviation of the performance and the mean computation time are
 * reported. 
 *
 * Finally, the overall best tasks are written to the specified output. These
 * tasks are selected to have both the highest mean performance, as well as the
 * smallest standard deviation in their performance. This is done as follows.
 * First the 99th percentile of task performance and the 1st percentile of
 * standard deviation is calculated. If a task exists for which the mean
 * performance of the repeats and the standard deviation equals these values
 * respectively, this task is found to be the best and is written to the
 * output. If no such task exists, the 98th percentile of performance and the
 * 2nd percentile of standard deviation is considered. This is repeated until
 * an interval is found which contains tasks. If one or more tasks are found,
 * this loop stops.
 *
 * @param[in] 	q 		Queue of Task structs which have already been 
 * 				run and have a Task::performance value
 * @param[in] 	repeats 	Number of times to repeat the best
 * 				configurations for consistency
 * @param[in] 	traintype 	type of training to do (CV or TT)
 *
 */
void consistency_repeats(struct Queue *q, long repeats, TrainType traintype)
{
	long i, f, r, N, *cv_idx;
	double p, pi, pr, pt, *time, *std, *mean, *perf;
	struct Queue *nq;
	struct GenData **train_folds, **test_folds;
	struct GenModel *model = gensvm_init_model();
	struct Task *task = NULL;
	clock_t loop_s, loop_e;

	nq = create_top_queue(q);
	N = nq->N;

	note("Number of items: %li\n", nq->N);
	std = Calloc(double, N);
	mean = Calloc(double, N);
	time = Calloc(double, N);
	perf = Calloc(double, N*repeats);

	task = get_next_task(nq);

	model->n = 0;
	model->m = task->train_data->m;
	model->K = task->train_data->K;
	gensvm_allocate_model(model);
	gensvm_seed_model_V(NULL, model, task->train_data);

	cv_idx = Calloc(long, task->train_data->n);

	train_folds = Malloc(struct GenData *, task->folds);
	test_folds = Malloc(struct GenData *, task->folds);

	i = 0;
	while (task) {
		make_model_from_task(task, model);

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

			loop_s = clock();
			p = gensvm_cross_validation(model, train_folds, test_folds,
					task->folds, task->train_data->n);
			loop_e = clock();
			time[i] += elapsed_time(loop_s, loop_e);
			matrix_set(perf, repeats, i, r, p);
			mean[i] += p/((double) repeats);
			note("%3.3f\t", p);
			// this is done because if we reuse the V it's not a 
			// consistency check
			gensvm_seed_model_V(NULL, model, task->train_data);
			for (f=0; f<task->folds; f++) {
				gensvm_free_data(train_folds[f]);
				gensvm_free_data(test_folds[f]);
			}
			free(train_folds);
			free(test_folds);
		}
		for (r=0; r<repeats; r++) {
			std[i] += pow(matrix_get(perf, repeats, i, r) - mean[i], 2.0);
		}
		if (r > 1) {
			std[i] /= ((double) repeats) - 1.0;
			std[i] = sqrt(std[i]);
		} else {
			std[i] = 0.0;
		}
		note("(m = %3.3f, s = %3.3f, t = %3.3f)\n", mean[i], std[i], time[i]);
		task = get_next_task(nq);
		i++;
	}

	// find the best overall configurations: those with high average
	// performance and low deviation in the performance
	note("\nBest overall configuration(s):\n");
	note("ID\tweights\tepsilon\t\tp\t\tkappa\t\tlambda\t\t"
			"mean_perf\tstd_perf\ttime_perf\n");
	p = 0.0;
	bool breakout = false;
	while (breakout == false) {
		pi = prctile(mean, N, (100.0-p));
		pr = prctile(std, N, p);
		pt = prctile(time, N, p);
		for (i=0; i<N; i++) 
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
			}
		p += 1.0;
	}

	free(cv_idx);
	// make sure no double free occurs with the copied kernelparam
	model->kernelparam = NULL;
	gensvm_free_model(model);
	free(nq->tasks);
	free(nq);
	free(perf);
	free(std);
	free(mean);
	free(time);
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
bool kernel_changed(struct Task *newtask, struct Task *oldtask)
{
	int i;
	if (oldtask == NULL)
		return true;
	if (newtask->kerneltype != oldtask->kerneltype) {
		return true;
	} else if (newtask->kerneltype == K_POLY) {
		for (i=0; i<3; i++)
			if (newtask->kernelparam[i] != oldtask->kernelparam[i])
				return true;
		return false;
	} else if (newtask->kerneltype == K_RBF) {
		if (newtask->kernelparam[0] != oldtask->kernelparam[0])
			return true;
		return false;
	} else if (newtask->kerneltype == K_SIGMOID) {
		for (i=0; i<2; i++)
			if (newtask->kernelparam[i] != oldtask->kernelparam[i])
				return true;
		return false;
	}
	return false;
}

/**
 * @brief Run the grid search for a Queue
 *
 * @details
 * Given a Queue of Task struct to be trained, a grid search is launched to
 * find the optimal parameter configuration. As is also done within
 * cross_validation(), the optimal weights of one parameter set are used as
 * initial estimates for GenModel::V in the next parameter set. Note that to
 * optimally exploit this feature of the optimization algorithm, the order in
 * which tasks are considered is important. This is considered in 
 * make_queue().
 * 
 * The performance found by cross validation is stored in the Task struct. 
 *
 * @param[in,out] 	q 	Queue with Task instances to run
 */

void start_training(struct Queue *q)
{
	int f, folds;
	double perf, current_max = 0;
	struct Task *task = get_next_task(q);
	struct Task *prevtask = NULL;
	struct GenModel *model = gensvm_init_model();
	clock_t main_s, main_e, loop_s, loop_e;

	// in principle this can change between tasks, but this shouldn't be 
	// the case TODO
	folds = task->folds;

	model->n = 0;
	model->m = task->train_data->m;
	model->K = task->train_data->K;
	gensvm_allocate_model(model);
	gensvm_seed_model_V(NULL, model, task->train_data);

	long *cv_idx = Calloc(long, task->train_data->n);
	gensvm_make_cv_split(task->train_data->n, task->folds, cv_idx);

	struct GenData **train_folds = Malloc(struct GenData *, task->folds);
	struct GenData **test_folds = Malloc(struct GenData *, task->folds);
	for (f=0; f<folds; f++) {
		train_folds[f] = gensvm_init_data();
		test_folds[f] = gensvm_init_data();
		gensvm_get_tt_split(task->train_data, train_folds[f],
			       	test_folds[f], cv_idx, f);
	}

	main_s = clock();
	while (task) {
		print_progress_string(task, q->N);
		make_model_from_task(task, model);
	
		if (kernel_changed(task, prevtask)) {
			note("*");
			for (f=0; f<folds; f++) {
				if (train_folds[f]->Z != train_folds[f]->RAW)
					free(train_folds[f]->Z);
				if (test_folds[f]->Z != test_folds[f]->RAW)
					free(test_folds[f]->Z);
				gensvm_kernel_preprocess(model,
					       	train_folds[f]);
				gensvm_kernel_postprocess(model,
					       	train_folds[f], test_folds[f]);
			}
			note("*");
		}

		loop_s = clock();
		perf = gensvm_cross_validation(model, train_folds, test_folds,
				folds, task->train_data->n);
		loop_e = clock();
		current_max = maximum(current_max, perf);

		note("\t%3.3f%% (%3.3fs)\t(best = %3.3f%%)\n", perf,
				elapsed_time(loop_s, loop_e), current_max);

		q->tasks[task->ID]->performance = perf;
		prevtask = task;
		task = get_next_task(q);
	}
	main_e = clock();

	note("\nTotal elapsed training time: %8.8f seconds\n",
			elapsed_time(main_s, main_e));

	// make sure no double free occurs with the copied kernelparam
	model->kernelparam = NULL;
	gensvm_free_model(model);
	for (f=0; f<folds; f++) {
		gensvm_free_data(train_folds[f]);
		gensvm_free_data(test_folds[f]);
	}
	free(train_folds);
	free(test_folds);
	free(cv_idx);
}

/**
 * @brief Run cross validation with a given set of train/test folds
 *
 * @details
 * This cross validation function uses predefined train/test splits. Also, the  
 * the optimal parameters GenModel::V of a previous fold as initial conditions 
 * for GenModel::V of the next fold. 
 *
 * @param[in] 	model 		GenModel with the configuration to train
 * @param[in] 	train_folds 	array of training datasets
 * @param[in] 	test_folds 	array of test datasets
 * @param[in] 	folds 		number of folds
 * @param[in] 	n_total 	number of objects in the union of the train 
 * datasets
 * @return 			performance (hitrate) of the configuration on 
 * cross validation
 */
double gensvm_cross_validation(struct GenModel *model,
	      	struct GenData **train_folds, struct GenData **test_folds,
		int folds, long n_total)
{
	FILE *fid;

	int f;
	long *predy;
	double performance, total_perf = 0;

	for (f=0; f<folds; f++) {
		// reallocate model in case dimensions differ with data
		gensvm_reallocate_model(model, train_folds[f]->n,
				train_folds[f]->r);

		// initialize object weights
		gensvm_initialize_weights(train_folds[f], model);

		// train the model (surpressing output)
		fid = GENSVM_OUTPUT_FILE;
		GENSVM_OUTPUT_FILE = NULL;
		gensvm_optimize(model, train_folds[f]);
		GENSVM_OUTPUT_FILE = fid;

		// calculate prediction performance on test set
		predy = Calloc(long, test_folds[f]->n);
		gensvm_predict_labels(test_folds[f], model, predy);
		performance = gensvm_prediction_perf(test_folds[f], predy);
		total_perf += performance * test_folds[f]->n;

		free(predy);
	}

	total_perf /= ((double) n_total);

	return total_perf;
}	


/**
 * @brief Free the Queue struct
 *
 * @details
 * Freeing the allocated memory of the Queue means freeing every Task struct
 * and then freeing the Queue.
 *
 * @param[in] 	q 	Queue to be freed
 *
 */
void free_queue(struct Queue *q)
{
	long i;
	for (i=0; i<q->N; i++) {
		free(q->tasks[i]->kernelparam);
		free(q->tasks[i]);
	}
	free(q->tasks);
	free(q);
}

/**
 * @brief Copy parameters from Task to GenModel
 *
 * @details
 * A Task struct only contains the parameters of the GenModel to be estimated.
 * This function is used to copy these parameters.
 *
 * @param[in] 		task 	Task instance with parameters
 * @param[in,out] 	model 	GenModel to which the parameters are copied
 */
void make_model_from_task(struct Task *task, struct GenModel *model)
{
	// copy basic model parameters
	model->weight_idx = task->weight_idx;
	model->epsilon = task->epsilon;
	model->p = task->p;
	model->kappa = task->kappa;
	model->lambda = task->lambda;

	// copy kernel parameters
	model->kerneltype = task->kerneltype;
	model->kernelparam = task->kernelparam;
}

/**
 * @brief Copy model parameters between two GenModel structs
 *
 * @details
 * The parameters copied are GenModel::weight_idx, GenModel::epsilon,
 * GenModel::p, GenModel::kappa, and GenModel::lambda.
 *
 * @param[in] 		from 	GenModel to copy parameters from
 * @param[in,out] 	to 	GenModel to copy parameters to
 */
void copy_model(struct GenModel *from, struct GenModel *to)
{
	to->weight_idx = from->weight_idx;
	to->epsilon = from->epsilon;
	to->p = from->p;
	to->kappa = from->kappa;
	to->lambda = from->lambda;

	to->kerneltype = from->kerneltype;
	switch (to->kerneltype) {
		case K_LINEAR:
			break;
		case K_POLY:
			to->kernelparam = Malloc(double, 3);
			to->kernelparam[0] = from->kernelparam[0];
			to->kernelparam[1] = from->kernelparam[1];
			to->kernelparam[2] = from->kernelparam[2];
			break;
		case K_RBF:
			to->kernelparam = Malloc(double, 1);
			to->kernelparam[0] = from->kernelparam[0];
			break;
		case K_SIGMOID:
			to->kernelparam = Malloc(double, 2);
			to->kernelparam[0] = from->kernelparam[0];
			to->kernelparam[1] = from->kernelparam[1];
			break;
	}
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
 * @param[in] 	task 	the Task specified
 * @param[in] 	N 	total number of tasks
 *
 */
void print_progress_string(struct Task *task, long N)
{
	char buffer[MAX_LINE_LENGTH];
	sprintf(buffer, "(%03li/%03li)\t", task->ID+1, N);
	if (task->kerneltype == K_POLY)
		sprintf(buffer + strlen(buffer), "d = %2.2f\t",
				task->kernelparam[2]);
	if (task->kerneltype == K_POLY || task->kerneltype == K_SIGMOID)
		sprintf(buffer + strlen(buffer), "c = %2.2f\t",
				task->kernelparam[1]);
	if (task->kerneltype == K_POLY || task->kerneltype == K_SIGMOID ||
			task->kerneltype == K_RBF)
		sprintf(buffer + strlen(buffer), "g = %3.3f\t",
				task->kernelparam[0]);
	sprintf(buffer + strlen(buffer), "eps = %g\tw = %i\tk = %2.2f\t"
			"l = %f\tp = %2.2f\t", task->epsilon,
			task->weight_idx, task->kappa, task->lambda, task->p);
	note(buffer);
}
