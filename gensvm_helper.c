#include <stdlib.h>
#include <numpy/arrayobject.h>
#include "include/gensvm_base.h"
#include "include/gensvm_predict.h"
#include "include/gensvm_task.h"
#include "include/gensvm_gridsearch.h"

void set_model(struct GenModel *model, double p, double lambda, 
		double kappa, double epsilon, int weight_idx, int kernel_index,
		double degree, double gamma, double coef, 
		double kernel_eigen_cutoff, long max_iter, long random_seed)
{
	model->p = p;
	model->lambda = lambda;
	model->kappa = kappa;
	model->epsilon = epsilon;
	model->weight_idx = weight_idx;
	model->kerneltype = kernel_index;
	model->gamma = gamma;
	model->coef = coef;
	model->degree = degree;
	model->kernel_eigen_cutoff = kernel_eigen_cutoff;
	model->max_iter = max_iter;
	model->seed = random_seed;
}


void set_seed_model(struct GenModel *model, double p, double lambda, 
		double kappa, double epsilon, int weight_idx, int kernel_index,
		double degree, double gamma, double coef, 
		double kernel_eigen_cutoff, long max_iter, long random_seed, 
		char *seed_V, long n_var, long n_class)
{
	model->p = p;
	model->lambda = lambda;
	model->kappa = kappa;
	model->epsilon = epsilon;
	model->weight_idx = weight_idx;
	model->kerneltype = kernel_index;
	model->gamma = gamma;
	model->coef = coef;
	model->degree = degree;
	model->kernel_eigen_cutoff = kernel_eigen_cutoff;
	model->max_iter = max_iter;
	model->seed = random_seed;

	model->n = 0;
	model->m = n_var;
	model->K = n_class;

	gensvm_allocate_model(model);

	double *seed_Vd = (double *) seed_V;

	long i, j;
	for (i=0; i<model->m+1; i++) {
		for (j=0; j<model->K-1; j++) {
			model->V[i*(n_class-1)+j] = seed_Vd[i*(n_class-1)+j];
		}
	}
}


void copy_X(struct GenData *data, double *Xd)
{
	// unfortunately, because we need the column of ones, we'll have to 
	// copy.
	long i, j;
	data->RAW = malloc(data->n * (data->m + 1) * sizeof(double));
	for (i=0; i<data->n; i++) {
		data->RAW[i*(data->m+1)] = 1.0;
		for (j=0; j<data->m; j++) {
			data->RAW[i*(data->m+1) + j+1] = Xd[i*data->m + j];
		}
	}
	data->Z = data->RAW;
}

void set_data(struct GenData *data, char *X, char *Y, npy_intp *dims, 
		long nr_class)
{
	if (data == NULL)
		return;
	data->n = (long) dims[0];
	data->m = (long) dims[1];
	data->r = data->m;
	data->K = nr_class;

	double *Xd = (double *) X;
	long *Yd = (long *) Y;

	copy_X(data, Xd);
	data->y = Yd;
}

const char *check_model(struct GenModel *model)
{
	if (model->epsilon <= 0)
		return "epsilon <= 0";
	if (model->kappa <= -1.0)
		return "kappa <= -1.0";
	if (model->lambda <= 0)
		return "lambda <= 0";
	if (model->p < 1.0 || model->p > 2.0)
		return "p not in [1, 2]";

	return NULL;
}

void copy_V(void *data, struct GenModel *model)
{
	memcpy(data, model->V, (model->m+1)*(model->K-1)*sizeof(double));
}

void copy_V_to_model(void *data, struct GenModel *model)
{
	memcpy(model->V, data, (model->m+1)*(model->K-1)*sizeof(double));
}

long get_iter_count(struct GenModel *model)
{
	return model->elapsed_iter;
}

double get_training_error(struct GenModel *model)
{
	return model->training_error;
}

int get_status(struct GenModel *model)
{
	return model->status;
}

long get_n(struct GenModel *model)
{
	return model->n;
}

long get_m(struct GenModel *model)
{
	return model->m;
}

long get_K(struct GenModel *model)
{
	return model->K;
}

// we don't free y, or data itself, because those were allocated by Python. 
void free_data(struct GenData *data)
{
	if (data == NULL)
		return;
	if (data->spZ != NULL)
		gensvm_free_sparse(data->spZ);
	if (data->Z == data->RAW) {
		free(data->Z);
	} else {
		free(data->Z);
		free(data->RAW);
	}
	free(data->Sigma);
}

void set_verbosity(int verbosity_flag)
{
	extern FILE *GENSVM_OUTPUT_FILE;
	extern FILE *GENSVM_ERROR_FILE;
	if (verbosity_flag) {
		GENSVM_OUTPUT_FILE = stdout;
		GENSVM_ERROR_FILE = stderr;
	} else {
		GENSVM_OUTPUT_FILE = NULL;
		GENSVM_ERROR_FILE = NULL;
	}
}

void gensvm_predict(char *X, char *V, long n_test, long m, long K, 
		char *predictions)
{
	struct GenModel *model = gensvm_init_model();
	struct GenData *data = gensvm_init_data();
	double *Xd = (double *) X;
	double *Vd = (double *) V;
	long *pred = (long *) predictions;

	data->n = n_test;
	data->m = m;
	data->r = m;
	data->K = K;

	copy_X(data, Xd);

	model->n = 1;
	model->m = m;
	model->K = K;

	gensvm_allocate_model(model);
	copy_V_to_model(Vd, model);

	gensvm_predict_labels(data, model, pred);

	gensvm_free_model(model);
	gensvm_free_data(data);
}

void set_task(struct GenTask *t, int ID, struct GenData *data, int folds, 
		double p, double lambda, double kappa, double epsilon, 
		double weight_idx, int kernel_index, double degree, 
		double gamma, double coef, long max_iter)
{
	t->ID = ID;
	t->train_data = data;
	t->folds = folds;
	t->p = p;
	t->lambda = lambda;
	t->kappa = kappa;
	t->epsilon = epsilon;
	t->weight_idx = weight_idx;
	t->kerneltype = kernel_index;
	t->degree = degree;
	t->gamma = gamma;
	t->coef = coef;
	t->max_iter = max_iter;
}

void gensvm_train_q_helper(struct GenQueue *q, char *CV_IDX, 
		int store_pred)
{
	long *cv_idx = (long *) CV_IDX;
	bool store_predictions = store_pred;

	gensvm_train_queue(q, cv_idx, store_predictions);
}

void set_queue(struct GenQueue *q, long n_tasks, struct GenTask **tasks)
{
	q->N = n_tasks;
	q->tasks = tasks;
}

double get_task_duration(struct GenTask *t)
{
	return t->duration;
}

double get_task_performance(struct GenTask *t)
{
	return t->performance;
}

void copy_task_predictions(struct GenTask *t, char *predictions, long n_obs)
{
	long *pred = (long *) predictions;
	long i;
	for (i=0; i<n_obs; i++) {
		pred[i] = t->predictions[i];
	}
}
