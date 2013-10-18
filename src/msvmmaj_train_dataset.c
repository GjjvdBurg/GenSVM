#include <math.h>
#include <time.h>

#include "crossval.h"
#include "libMSVMMaj.h"
#include "matrix.h"
#include "msvmmaj_train.h"
#include "msvmmaj_train_dataset.h"
#include "msvmmaj_pred.h"
#include "MSVMMaj.h"
#include "util.h"
#include "timer.h"

extern FILE *MSVMMAJ_OUTPUT_FILE;

void make_queue(struct Training *training, struct Queue *queue, 
		struct MajData *train_data, struct MajData *test_data)
{
	long i, j, k, l, m;
	long N, cnt = 0;
	struct Task *task;
	queue->i = 0;

	N = training->Np;
	N *= training->Nl;
	N *= training->Nk;
	N *= training->Ne;
	N *= training->Nw;

	queue->tasks = Malloc(struct Task *, N);
	queue->N = N;

	for (i=0; i<training->Ne; i++)
		for (j=0; j<training->Nw; j++)
			for (k=0; k<training->Nk; k++)
				for (l=0; l<training->Nl; l++)
					for (m=0; m<training->Np; m++) {
						task = Malloc(struct Task, 1);
						task->epsilon = training->epsilons[i];
						task->weight_idx = training->weight_idxs[j];
						task->kappa = training->kappas[k];
						task->lambda = training->lambdas[l];
						task->p = training->ps[m];
						task->train_data = train_data;
						task->test_data = test_data;
						task->folds = training->folds;
						task->ID = cnt;
						queue->tasks[cnt] = task;
						cnt++;
					}
}

struct Task *get_next_task(struct Queue *q)
{
	long i = q->i;
	if (i < q->N) {
		q->i++;
		return q->tasks[i];
	}
	return NULL;
}

int tasksort(const void *elem1, const void *elem2)
{
	const struct Task *t1 = (*(struct Task **) elem1);
	const struct Task *t2 = (*(struct Task **) elem2);
	return (t1->performance > t2->performance);
}

int doublesort(const void *elem1, const void *elem2)
{
	const double t1 = (*(double *) elem1);
	const double t2 = (*(double *) elem2);
	return t1 > t2;
}


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

void consistency_repeats(struct Queue *q, long repeats, TrainType traintype)
{
	long i, r, N;
	double p, pi, pr, boundary, time, *std, *mean, *perf;
	struct Queue *nq = Malloc(struct Queue, 1);
	struct MajModel *model = Malloc(struct MajModel, 1);
	struct Task *task = Malloc(struct Task, 1);
	clock_t loop_s, loop_e;

	// calculate the percentile (Matlab style)
	qsort(q->tasks, q->N, sizeof(struct Task *), tasksort);
	p = 0.95*q->N + 0.5;
	pi = maximum(minimum(floor(p), q->N-1), 1);
	pr = maximum(minimum(p - pi, 1), 0);
	boundary = (1 - pr)*q->tasks[((long) pi)-1]->performance;
	boundary += pr*q->tasks[((long) pi)]->performance;
	note("boundary determined at: %f\n", boundary);

	N = 0;
	for (i=0; i<q->N; i++)
		if (q->tasks[i]->performance >= boundary)
			N++;
	note("Number of items: %li\n", N);
	std = Calloc(double, N);
	mean = Calloc(double, N);
	perf = Calloc(double, N*repeats);

	nq->tasks = Malloc(struct Task *, N);
	for (i=q->N-1; i>q->N-N-1; i--) 
		nq->tasks[q->N-i-1] = q->tasks[i];
	nq->N = N;
	nq->i = 0;

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
				p = cross_validation(model, NULL, task->train_data, task->folds);
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
			std[i] += pow(matrix_get(perf, repeats, i, r) - mean[i], 2);
		}
		std[i] /= ((double) repeats) - 1.0;
		std[i] = sqrt(std[i]);
		note("(m = %3.3f, s = %3.3f, t = %3.3f)\n", mean[i], std[i], time);
	}

	note("\nBest overall configuration(s):\n");
	note("ID\tweights\tepsilon\t\tp\t\tkappa\t\tlambda\t\tmean_perf\tstd_perf\n");
	p = 0.0;
	bool breakout = false;
	while (breakout == false) {
		pi = prctile(mean, N, (100.0-p)/100.0);
		pr = prctile(std, N, p/100.0);
		for (i=0; i<N; i++) 
			if ((pi - mean[i] < 0.0001) && (std[i] - pr < 0.0001)) {
				note("(%li)\tw = %li\te = %f\tp = %f\tk = %f\tl = %f\t"
						"mean: %3.3f\tstd: %3.3f\n",
						nq->tasks[i]->ID, 
						nq->tasks[i]->weight_idx,
						nq->tasks[i]->epsilon, nq->tasks[i]->p,
						nq->tasks[i]->kappa, nq->tasks[i]->lambda,
						mean[i], std[i]);
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
		seed_model = Malloc(struct MajModel, 1);
		seed_model->n = 0; // we never use anything other than V
		seed_model->m = model->m;
		seed_model->K = model->K;
		msvmmaj_allocate_model(seed_model);
		msvmmaj_seed_model_V(NULL, seed_model);
		fs = true;
	}

	train_data = Malloc(struct MajData, 1);
	test_data = Malloc(struct MajData, 1);

	msvmmaj_make_cv_split(model->n, folds, cv_idx);
	for (f=0; f<folds; f++) {
		msvmmaj_get_tt_split(data, train_data, test_data, cv_idx, f);

		fold_model = Malloc(struct MajModel, 1);
		copy_model(model, fold_model);
	
		fold_model->n = train_data->n;
		fold_model->m = train_data->m;
		fold_model->K = train_data->K;
		
		msvmmaj_allocate_model(fold_model);
		msvmmaj_initialize_weights(train_data, fold_model);
		msvmmaj_seed_model_V(seed_model, fold_model);
	
		fid = MSVMMAJ_OUTPUT_FILE;
		MSVMMAJ_OUTPUT_FILE = NULL;
		msvmmaj_optimize(fold_model, train_data);
		MSVMMAJ_OUTPUT_FILE = fid;

		predy = Calloc(long, test_data->n);
		msvmmaj_predict_labels(test_data, fold_model, predy);
		performance[f] = msvmmaj_prediction_perf(test_data, predy);
		total_perf += performance[f]/((double) folds);

		msvmmaj_seed_model_V(fold_model, seed_model);
		
		free(predy);
		free(train_data->y);
		free(train_data->Z);
		free(test_data->y);
		free(test_data->Z);

		msvmmaj_free_model(fold_model);
	}

	if (fs)
		msvmmaj_free_model(seed_model);
	free(train_data);
	free(test_data);
	free(performance);
	free(cv_idx);

	return total_perf;

}

void start_training_cv(struct Queue *q)
{
	double perf, current_max = 0;
	struct Task *task = get_next_task(q);
	struct MajModel *seed_model = Malloc(struct MajModel, 1);
	struct MajModel *model = Malloc(struct MajModel, 1);
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
		note("(%03li/%03li)\tw = %li\te = %f\tp = %f\tk = %f\t l = %f\t",
				task->ID+1, q->N, task->weight_idx, task->epsilon,
				task->p, task->kappa, task->lambda);
		make_model_from_task(task, model);
		
		loop_s = clock();
		perf = cross_validation(model, seed_model, task->train_data, task->folds);
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

void start_training_tt(struct Queue *q)
{
	FILE *fid;

	long c = 0;
	long *predy;
	double total_perf, current_max = 0;

	struct Task *task = get_next_task(q);
	struct MajModel *seed_model = Malloc(struct MajModel, 1);

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
		struct MajModel *model = Malloc(struct MajModel, 1);
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

void free_queue(struct Queue *q)
{
	long i;
	for (i=0; i<q->N; i++)
		free(q->tasks[i]);
	free(q->tasks);
	free(q);
}

void make_model_from_task(struct Task *task, struct MajModel *model)
{
	model->weight_idx = task->weight_idx;
	model->epsilon = task->epsilon;
	model->p = task->p;
	model->kappa = task->kappa;
	model->lambda = task->lambda;
}

void copy_model(struct MajModel *from, struct MajModel *to)
{
	to->weight_idx = from->weight_idx;
	to->epsilon = from->epsilon;
	to->p = from->p;
	to->kappa = from->kappa;
	to->lambda = from->lambda;
}
