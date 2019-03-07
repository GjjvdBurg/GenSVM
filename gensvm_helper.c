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

void set_raw_weights(struct GenModel *model, char *raw_weights, int n_obs)
{
	double *weights = (double *)raw_weights;
	if (model->rho != NULL)
		free(model->rho);
	model->rho = Calloc(double, n_obs);
	int i;
	for (i=0; i<n_obs; i++)
		model->rho[i] = weights[i];
}

struct GenData *_build_gensvm_data(double *X, int *y, int n, int m, int K)
{
	int i, j;
	double value;

	struct GenData *data = gensvm_init_data();
	data->n = n;
	data->m = m;
	data->r = m;
	data->K = K;

	data->RAW = Calloc(double, n*(m+1));

	for (i=0; i<n; i++) {
		for (j=0; j<m; j++) {
			value = matrix_get(X, n, m, i, j);
			matrix_set(data->RAW, n, m+1, i, j+1, value);
		}
		matrix_set(data->RAW, n, m+1, i, 0, 1.0);
	}
	data->Z = data->RAW;

	// convert to sparse matrix if possible
	if (gensvm_could_sparse(data->Z, n, m+1)) {
		note("Converting to sparse ... ");
		data->spZ = gensvm_dense_to_sparse(data->Z, n, m+1);
		note("done.\n");
		free(data->RAW);
		data->RAW = NULL;
		data->Z = NULL;
	}

	if (y == NULL) {
		data->y = NULL;
	} else {
		data->y = Malloc(long, n);
		for (i=0; i<n; i++)
			data->y[i] = y[i];
	}

	return data;
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

// we don't free y, because that was allocated by Python. 
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
	free(data);
}

void set_verbosity(int verbosity)
{
	extern FILE *GENSVM_OUTPUT_FILE;
	extern FILE *GENSVM_ERROR_FILE;
	extern void (*gensvm_print_out)(const char *, ...);
	extern void (*gensvm_print_err)(const char *, ...);

	gensvm_print_out = gensvm_print_output_fpt;
	gensvm_print_err = gensvm_print_error_fpt;

	if (verbosity) {
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
	double *Xd = (double *) X;
	double *Vd = (double *) V;
	long *pred = (long *) predictions;

	struct GenModel *model = gensvm_init_model();
	model->n = 1;
	model->m = m;
	model->K = K;
	gensvm_allocate_model(model);
	copy_V_to_model(Vd, model);

	struct GenData *data = _build_gensvm_data(Xd, NULL, n_test, m, K);

	gensvm_predict_labels(data, model, pred);

	gensvm_free_model(model);
	gensvm_free_data(data);
}

void gensvm_predict_kernels(char *X_test, char *X_train, char *V, long V_row, 
		long V_col, long n_train, long n_test, long m, long K, 
		int kernel_idx, double gamma, double coef, double degree,
		double kernel_eigen_cutoff, char *predictions)
{
	long i, j;
	double value;

	struct GenModel *model = gensvm_init_model();

	double *Xtrain = (double *) X_train;
	double *Xtest = (double *) X_test;
	double *Vd = (double *) V;
	long *pred = (long *) predictions;

	model->n = n_train;
	model->m = V_row - 1;
	model->K = V_col + 1;
	model->kerneltype = kernel_idx;
	model->gamma = gamma;
	model->coef = coef;
	model->degree = degree;
	model->kernel_eigen_cutoff = kernel_eigen_cutoff;

	gensvm_allocate_model(model);

	struct GenData *traindata = _build_gensvm_data(Xtrain, NULL, n_train,
			m, K);
	struct GenData *testdata = _build_gensvm_data(Xtest, NULL, n_test, m,
			K);

	gensvm_kernel_preprocess(model, traindata);
	gensvm_reallocate_model(model, traindata->n, traindata->r);

	for (i=0; i<model->m+1; i++) {
		for (j=0; j<model->K; j++) {
			value = matrix_get(Vd, V_row, V_col, i, j);
			matrix_set(model->V, model->m+1, model->K-1, i, j, value);
		}
	}

	gensvm_kernel_postprocess(model, traindata, testdata);

	gensvm_predict_labels(testdata, model, pred);

	gensvm_free_data(traindata);
	gensvm_free_data(testdata);
	gensvm_free_model(model);
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
		int store_pred, int verbosity)
{
	long *cv_idx = (long *) CV_IDX;
	bool store_predictions = store_pred;

	gensvm_train_queue(q, cv_idx, store_predictions, verbosity);
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
	long i, val;
	if (t->predictions == NULL) { // interrupt occured
		for (i=0; i<n_obs; i++)
			pred[i] = -1;
	} else {
		for (i=0; i<n_obs; i++) {
			val = t->predictions[i];
			pred[i] = val;
		}
	}
}

void copy_task_durations(struct GenTask *t, char *durations, int n_folds)
{
	double val, *dur = (double *)durations;
	int i;
	if (t->durations == NULL) { // interrupt occured
		for (i=0; i<n_folds; i++)
			dur[i] = NPY_NAN;
	} else {
		for (i=0; i<n_folds; i++) {
			val = t->durations[i];
			dur[i] = (val < 0) ? NPY_NAN : val;
		}
	}
}
