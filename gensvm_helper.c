#include <stdlib.h>
#include <numpy/arrayobject.h>
#include "include/gensvm_base.h"
#include "include/gensvm_predict.h"

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
