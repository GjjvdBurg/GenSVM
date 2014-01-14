/**
 * @file msvmmaj_train_dataset.c
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Functions for finding the optimal parameters for the dataset
 *
 * @details
 * The MSVMMaj algorithm takes a number of parameters. The functions in
 * this file are used to find the optimal parameters. 
 */

#include <math.h>
#include <time.h>

#include "crossval.h"
#include "libMSVMMaj.h"
#include "msvmmaj.h"
#include "msvmmaj_init.h"
#include "msvmmaj_matrix.h"
#include "msvmmaj_train.h"
#include "msvmmaj_train_dataset.h"
#include "msvmmaj_pred.h"
#include "util.h"
#include "timer.h"

extern FILE *MSVMMAJ_OUTPUT_FILE;

/**
 * @brief Initialize a Queue from a Training instance
 *
 * @details
 * A Training instance describes the grid to search over. This funtion
 * creates all tasks that need to be performed and adds these to 
 * a Queue. Each task contains a pointer to the train and test datasets
 * which are supplied. Note that the tasks are created in a specific order of
 * the parameters, to ensure that the MajModel::V of a previous parameter
 * set provides the best possible initial estimate of MajModel::V for the next
 * parameter set.
 *
 * @param[in] 	training 	Training struct describing the grid search
 * @param[in] 	queue 		pointer to a Queue that will be used to 
 * 				add the tasks to
 * @param[in] 	train_data 	MajData of the training set
 * @param[in] 	test_data 	MajData of the test set
 *
 */
void make_queue(struct Training *training, struct Queue *queue, 
		struct MajData *train_data, struct MajData *test_data)
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

	queue->tasks = Malloc(struct Task *, N);
	queue->N = N;

	// initialize all tasks
	for (i=0; i<N; i++) {
		task = Malloc(struct Task, 1);
		task->ID = i;
		task->train_data = train_data;
		task->test_data = test_data;
		task->folds = training->folds;
		task->kerneltype = training->kerneltype;
		task->kernel_param = Calloc(double, training->Ng + 
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
				queue->tasks[i]->kernel_param[0] = 
					training->gammas[j];
				i++;
			}

	cnt *= training->Ng > 0 ? training->Ng : 1;
	i = 0;
	while (i < N && training->Nc > 0)
		for (j=0; j<training->Nc; j++)
			for (k=0; k<cnt; k++) {
				queue->tasks[i]->kernel_param[1] = 
					training->coefs[j];
				i++;
			}

	cnt *= training->Nc > 0 ? training->Ng : 1;
	i = 0;
	while (i < N && training->Nd > 0)
		for (j=0; j<training->Nd; j++)
			for (k=0; k<cnt; k++) {
				queue->tasks[i]->kernel_param[2] = 
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
 * @brief Comparison function for doubles
 *
 * @details
 * Similar to tasksort() only now for two doubles.
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
 * @param[in] 	p 	percentile to calculate ( 0 <= p <= 1.0 ).
 * @returns 		the p-th percentile of the values
 */
double prctile(double *values, long N, double p)
{
	long i;
	double pi, pr, boundary;
	double *local = Malloc(double, N);
	for (i=0; i<N; i++)
		local[i] = values[i];

	qsort(local, N, sizeof(double), doublesort);
	p = p*N + 0.5;
	pi = maximum(minimum(floor(p), N-1), 1);
	pr = maximum(minimum(p - pi, 1), 0);
	boundary = (1 - pr)*local[((long) pi)-1] + pr*local[((long) pi)];

	free(local);

	return boundary;
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
	long i, r, N;
	double p, pi, pr, boundary, time, *std, *mean, *perf;
	struct Queue *nq = Malloc(struct Queue, 1);
	struct MajModel *model = msvmmaj_init_model();
	struct Task *task = Malloc(struct Task, 1);
	clock_t loop_s, loop_e;

	// calculate the performance percentile (Matlab style)
	qsort(q->tasks, q->N, sizeof(struct Task *), tasksort);
	p = 0.95*q->N + 0.5;
	pi = maximum(minimum(floor(p), q->N-1), 1);
	pr = maximum(minimum(p - pi, 1), 0);
	boundary = (1 - pr)*q->tasks[((long) pi)-1]->performance;
	boundary += pr*q->tasks[((long) pi)]->performance;
	note("boundary determined at: %f\n", boundary);
	
	// find the number of tasks that perform at least as good as the 95th
	// percentile
	N = 0;
	for (i=0; i<q->N; i++)
		if (q->tasks[i]->performance >= boundary)
			N++;
	note("Number of items: %li\n", N);
	std = Calloc(double, N);
	mean = Calloc(double, N);
	perf = Calloc(double, N*repeats);

	// create a new task queue with the tasks which perform well
	nq->tasks = Malloc(struct Task *, N);
	for (i=q->N-1; i>q->N-N-1; i--) 
		nq->tasks[q->N-i-1] = q->tasks[i];
	nq->N = N;
	nq->i = 0;

	// for each task run the consistency repeats
	for (i=0; i<N; i++) {
		task = get_next_task(nq);
		make_model_from_task(task, model);

		model->n = task->train_data->n;
		model->m = task->train_data->m;
		model->K = task->train_data->K;

		time = 0;
		note("(%02li/%02li:%03li)\t", i+1, N, task->ID);
		for (r=0; r<repeats; r++) {
			if (traintype == CV) {
				loop_s = clock();
				p = cross_validation(model, NULL, 
						task->train_data, task->folds);
				loop_e = clock();
				time += elapsed_time(loop_s, loop_e);
				matrix_set(perf, repeats, i, r, p);
				mean[i] += p/((double) repeats);
			} else {
				note("Only cv is implemented\n");
				exit(1);
			}
			note("%3.3f\t", p);
		}
		for (r=0; r<repeats; r++) {
			std[i] += pow(matrix_get(
						perf, 
						repeats, 
						i, 
						r) - mean[i],
				       	2.0);
		}
		std[i] /= ((double) repeats) - 1.0;
		std[i] = sqrt(std[i]);
		note("(m = %3.3f, s = %3.3f, t = %3.3f)\n", 
				mean[i], std[i], time);
	}

	// find the best overall configurations: those with high average
	// performance and low deviation in the performance
	note("\nBest overall configuration(s):\n");
	note("ID\tweights\tepsilon\t\tp\t\tkappa\t\tlambda\t\t"
			"mean_perf\tstd_perf\n");
	p = 0.0;
	bool breakout = false;
	while (breakout == false) {
		pi = prctile(mean, N, (100.0-p)/100.0);
		pr = prctile(std, N, p/100.0);
		for (i=0; i<N; i++) 
			if ((pi - mean[i] < 0.0001) && (std[i] - pr < 0.0001)) {
				note("(%li)\tw = %li\te = %f\tp = %f\t"
						"k = %f\tl = %f\t"
						"mean: %3.3f\tstd: %3.3f\n",
						nq->tasks[i]->ID, 
						nq->tasks[i]->weight_idx,
						nq->tasks[i]->epsilon, 
						nq->tasks[i]->p,
						nq->tasks[i]->kappa, 
						nq->tasks[i]->lambda,
						mean[i], 
						std[i]);
				breakout = true;
			}
		p += 1.0;
	}

	free(task);
	free(model);
	free(perf);
	free(std);
	free(mean);
}

/**
 * @brief Run cross validation with a seed model
 *
 * @details
 * This is an implementation of cross validation which uses the optimal
 * parameters MajModel::V of a previous fold as initial conditions for
 * MajModel::V of the next fold. An initial seed for V can be given through the
 * seed_model parameter. If seed_model is NULL, random starting values are
 * used.
 *
 * @todo
 * The seed model shouldn't have to be allocated completely, since only V is
 * used.
 * @todo
 * There must be some inefficiencies here because the fold model is allocated
 * at every fold. This would be detrimental with large datasets.
 *
 * @param[in] 	model 		MajModel with the configuration to train
 * @param[in] 	seed_model 	MajModel with a seed for MajModel::V
 * @param[in] 	data 		MajData with the dataset
 * @param[in] 	folds 		number of cross validation folds
 * @returns 			performance (hitrate) of the configuration on 
 * 				cross validation
 */
double cross_validation(struct MajModel *model, struct MajModel *seed_model,
		struct MajData *data, long folds) 
{
	FILE *fid;

	bool fs = false;
	long f, *predy;
	double total_perf = 0;
	struct MajModel *fold_model;
	struct MajData *train_data, *test_data;

	long *cv_idx = Calloc(long, model->n);
	double *performance = Calloc(double, folds);

	if (seed_model == NULL) {
		seed_model = msvmmaj_init_model();
		seed_model->n = 0; // we never use anything other than V
		seed_model->m = model->m;
		seed_model->K = model->K;
		msvmmaj_allocate_model(seed_model);
		msvmmaj_seed_model_V(NULL, seed_model);
		fs = true;
	}

	train_data = msvmmaj_init_data();
	test_data = msvmmaj_init_data();
 	// create splits
	msvmmaj_make_cv_split(model->n, folds, cv_idx);

	for (f=0; f<folds; f++) {
		msvmmaj_get_tt_split(data, train_data, test_data, cv_idx, f);
		// initialize a model for this fold and copy the model
		// parameters
		fold_model = msvmmaj_init_model();
		copy_model(model, fold_model);
	
		fold_model->n = train_data->n;
		fold_model->m = train_data->m;
		fold_model->K = train_data->K;
	
		// allocate, initialize and seed the fold model	
		msvmmaj_allocate_model(fold_model);
		msvmmaj_initialize_weights(train_data, fold_model);
		msvmmaj_seed_model_V(seed_model, fold_model);

		// train the model (without output)	
		fid = MSVMMAJ_OUTPUT_FILE;
		MSVMMAJ_OUTPUT_FILE = NULL;
		msvmmaj_optimize(fold_model, train_data);
		MSVMMAJ_OUTPUT_FILE = fid;

		// calculate predictive performance on test set 
		predy = Calloc(long, test_data->n);
		msvmmaj_predict_labels(test_data, fold_model, predy);
		performance[f] = msvmmaj_prediction_perf(test_data, predy);
		total_perf += performance[f]/((double) folds);

		// seed the seed model with the fold model
		msvmmaj_seed_model_V(fold_model, seed_model);
		
		free(predy);
		free(train_data->y);
		free(train_data->Z);
		free(test_data->y);
		free(test_data->Z);

		msvmmaj_free_model(fold_model);
	}

	// if a seed model was allocated before, free it.
	if (fs)
		msvmmaj_free_model(seed_model);
	free(train_data);
	free(test_data);
	free(performance);
	free(cv_idx);

	return total_perf;

}

/**
 * @brief Run the grid search for a cross validation dataset
 *
 * @details
 * Given a Queue of Task struct to be trained, a grid search is launched to
 * find the optimal parameter configuration. As is also done within
 * cross_validation(), the optimal weights of one parameter set are used as
 * initial estimates for MajModel::V in the next parameter set. Note that to
 * optimally exploit this feature of the optimization algorithm, the order in
 * which tasks are considered is important. This is considered in 
 * make_queue().
 * 
 * The performance found by cross validation is stored in the Task struct. 
 *
 * @param[in,out] 	q 	Queue with Task instances to run
 */
void start_training_cv(struct Queue *q)
{
	double perf, current_max = 0;
	struct Task *task = get_next_task(q);
	struct MajModel *seed_model = msvmmaj_init_model();
	struct MajModel *model = msvmmaj_init_model();
	clock_t main_s, main_e, loop_s, loop_e;

	model->n = task->train_data->n;
	model->m = task->train_data->m;
	model->K = task->train_data->K;
	msvmmaj_allocate_model(model);

	seed_model->n = 0;
	seed_model->m = task->train_data->m;
	seed_model->K = task->train_data->K;
	msvmmaj_allocate_model(seed_model);
	msvmmaj_seed_model_V(NULL, seed_model);
	
	main_s = clock();
	while (task) {
		note("(%03li/%03li)\tw = %li\te = %f\tp = %f\tk = %f\t "
				"l = %f\t",
				task->ID+1, q->N, task->weight_idx,
			       	task->epsilon,
				task->p, task->kappa, task->lambda);
		make_model_from_task(task, model);
		
		loop_s = clock();
		perf = cross_validation(model, seed_model, task->train_data,
			       	task->folds);
		loop_e = clock();
		current_max = maximum(current_max, perf);

		note("\t%3.3f%% (%3.3fs)\t(best = %3.3f%%)\n", perf,
				elapsed_time(loop_s, loop_e),
				current_max);

		q->tasks[task->ID]->performance = perf;
		task = get_next_task(q);
	}
	main_e = clock();

	note("\nTotal elapsed time: %8.8f seconds\n", 
			elapsed_time(main_s, main_e));

	free(task);
	msvmmaj_free_model(seed_model);
}

/**
 * @brief Run the grid search for a train/test dataset
 *
 * @details
 * This function is similar to start_training_cv(), except that the
 * pre-determined training set is used only once, and the pre-determined test
 * set is used for validation.
 *
 * @todo
 * It would probably be better to train the model on the training set using
 * cross validation and only use the test set when comparing with other
 * methods. The way it is now, you're finding out which parameters predict
 * _this_ test set best, which is not what you want.
 *
 * @param[in] 	q 	Queue with Task structs to run
 *
 */
void start_training_tt(struct Queue *q)
{
	FILE *fid;

	long c = 0;
	long *predy;
	double total_perf, current_max = 0;

	struct Task *task = get_next_task(q);
	struct MajModel *seed_model = msvmmaj_init_model();

	clock_t main_s, main_e;
	clock_t loop_s, loop_e;

	seed_model->m = task->train_data->m;
	seed_model->K = task->train_data->K;
	msvmmaj_allocate_model(seed_model);
	msvmmaj_seed_model_V(NULL, seed_model);
	
	main_s = clock();
	while (task) {
		total_perf = 0;
		note("(%li/%li)\tw = %li\te = %f\tp = %f\tk = %f\tl = %f\t",
				c+1, q->N, task->weight_idx, task->epsilon,
				task->p, task->kappa, task->lambda);
		loop_s = clock();
		struct MajModel *model = msvmmaj_init_model();
		make_model_from_task(task, model);

		model->n = task->train_data->n;
		model->m = task->train_data->m;
		model->K = task->train_data->K;

		msvmmaj_allocate_model(model);
		msvmmaj_initialize_weights(task->train_data, model);
		msvmmaj_seed_model_V(seed_model, model);

		fid = MSVMMAJ_OUTPUT_FILE;
		MSVMMAJ_OUTPUT_FILE = NULL;
		msvmmaj_optimize(model, task->train_data);
		MSVMMAJ_OUTPUT_FILE = fid;

		predy = Calloc(long, task->test_data->n);
		msvmmaj_predict_labels(task->test_data, model, predy);
		if (task->test_data->y != NULL) 
			total_perf = msvmmaj_prediction_perf(task->test_data, predy);
		msvmmaj_seed_model_V(model, seed_model);

		msvmmaj_free_model(model);
		free(predy);
		note(".");
		loop_e = clock();
		current_max = maximum(current_max, total_perf);
		note("\t%3.3f%% (%3.3fs)\t(best = %3.3f%%)\n", total_perf,
				elapsed_time(loop_s, loop_e), current_max);
		q->tasks[task->ID]->performance = total_perf;
		task = get_next_task(q);
	}
	main_e = clock();
	
	note("\nTotal elapsed time: %8.8f seconds\n",
			elapsed_time(main_s, main_e));
	free(task);
	msvmmaj_free_model(seed_model);
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
		free(q->tasks[i]->kernel_param);
		free(q->tasks[i]);
	}
	free(q->tasks);
	free(q);
}

/**
 * @brief Copy parameters from Task to MajModel
 *
 * @details
 * A Task struct only contains the parameters of the MajModel to be estimated.
 * This function is used to copy these parameters.
 *
 * @param[in] 		task 	Task instance with parameters
 * @param[in,out] 	model 	MajModel to which the parameters are copied
 */
void make_model_from_task(struct Task *task, struct MajModel *model)
{
	model->weight_idx = task->weight_idx;
	model->epsilon = task->epsilon;
	model->p = task->p;
	model->kappa = task->kappa;
	model->lambda = task->lambda;
}

/**
 * @brief Copy model parameters between two MajModel structs
 *
 * @details
 * The parameters copied are MajModel::weight_idx, MajModel::epsilon,
 * MajModel::p, MajModel::kappa, and MajModel::lambda.
 *
 * @param[in] 		from 	MajModel to copy parameters from
 * @param[in,out] 	to 	MajModel to copy parameters to
 */
void copy_model(struct MajModel *from, struct MajModel *to)
{
	to->weight_idx = from->weight_idx;
	to->epsilon = from->epsilon;
	to->p = from->p;
	to->kappa = from->kappa;
	to->lambda = from->lambda;
}
