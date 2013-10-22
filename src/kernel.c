/**
 * @file kernel.c
 * @author Gertjan van den Burg (burg@ese.eur.nl)
 * @date October 18, 2013
 * @brief Defines main functions for use of kernels in MSVMMaj.
 *
 * @details
 * Functions for constructing different kernels using user-supplied 
 * parameters. Also contains the functions for decomposing the
 * kernel matrix using several decomposition methods.
 *
 */
#include <math.h>

#include "kernel.h"

void msvmmaj_make_kernel(struct MajModel *model, struct MajData *data)
{
	switch (model->kerneltype) {
		case K_LINEAR:
			break;
		case K_POLY:
			msvmmaj_make_kernel_poly(model, data);
			break;
		case K_RBF:
			msvmmaj_make_kernel_rbf(model, data);
			break;
		case K_SIGMOID:
			msvmmaj_make_kernel_sigmoid(model, data);
			break;
	}
}

void msvmmaj_make_kernel_rbf(struct MajModel *model, struct MajData *data)
{
	long i, j;
	long n = model->n;
	double value;
	double *x1, *x2;
	double *K = Calloc(double, n*(n+1));

	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			x1 = &data->Z[i*(data->m+1)+1];
			x2 = &data->Z[j*(data->m+1)+1];
			value = msvmmaj_compute_rbf(x1, x2, model->kernelparam, n);
			matrix_set(K, n+1, i, j+1, value);
		}
		matrix_set(K, n+1, i, 0, 1.0);
	}

	free(data->Z);
	data->Z = K;
	data->m = n;
	model->m = n;
}

/**
 * Implements k(x, z) = exp( -gamma * || x - z ||^2)
 */
double msvmmaj_compute_rbf(double *x1, double *x2, double *kernelparam, long n)
{
	long i;
	double value = 0.0;

	for (i=0; i<n; i++) 
		value += (x1[i] - x2[i]) * (x1[i] - x2[i]);
	value *= -kernelparam[0];
	return exp(value);
}

/**
 * Implements k(x, z) = (gamma * <x, z> + c)^degree
 */
double msvmmaj_compute_poly(double *x1, double *x2, double *kernelparam, long n)
{
	long i;
	double value = 0.0;
	for (i=0; i<n; i++)
		value += x1[i]*x2[i];
	value *= kernelparam[0];
	value += kernelparam[1];
	for (i=1; i<(int kernelparam[2]); i++)
		value *= value;
	:w
