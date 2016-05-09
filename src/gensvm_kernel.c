/**
 * @file gensvm_kernel.c
 * @author Gertjan van den Burg
 * @date October 18, 2013
 * @brief Defines main functions for use of kernels in GenSVM.
 *
 * @details
 * Functions for constructing different kernels using user-supplied
 * parameters. Also contains the functions for decomposing the
 * kernel matrix using several decomposition methods.
 *
 */

#include <cblas.h>
#include <math.h>

#include "gensvm.h"
#include "gensvm_kernel.h"
#include "gensvm_lapack.h"
#include "gensvm_matrix.h"
#include "gensvm_util.h"

/**
 * @brief Do the preprocessing steps needed to perform kernel GenSVM
 *
 * @details
 * tdb
 *
 */
void gensvm_kernel_preprocess(struct GenModel *model, struct GenData *data)
{
	if (model->kerneltype == K_LINEAR) {
		data->r = data->m;
		return;
	}

	int i;
	long r,
	     n = data->n;
	double *P = NULL,
	       *Sigma = NULL,
	       *K = NULL;

	// build the kernel matrix
	K = Calloc(double, n*n);
	gensvm_make_kernel(model, data, K);

	// generate the eigen decomposition
	r = gensvm_make_eigen(K, n, &P, &Sigma);

	// build M and set to data (leave RAW intact)
	gensvm_make_trainfactor(data, P, Sigma, r);

	// Set Sigma to data->Sigma (need it again for prediction)
	if (data->Sigma != NULL)
		free(data->Sigma);
	data->Sigma = Sigma;

	// write kernel params to data
	data->kerneltype = model->kerneltype;
	free(data->kernelparam);
	switch (model->kerneltype) {
		case K_LINEAR:
			break;
		case K_POLY:
			data->kernelparam = Calloc(double, 3);
			for (i=0; i<3; i++)
				data->kernelparam[i] = model->kernelparam[i];
			break;
		case K_RBF:
			data->kernelparam = Calloc(double, 1);
			data->kernelparam[0] = model->kernelparam[0];
			break;
		case K_SIGMOID:
			data->kernelparam = Calloc(double, 2);
			data->kernelparam[0] = model->kernelparam[0];
			data->kernelparam[1] = model->kernelparam[1];
	}

	free(K);
	free(P);
}

void gensvm_kernel_postprocess(struct GenModel *model,
		struct GenData *traindata, struct GenData *testdata)
{
	if (model->kerneltype == K_LINEAR) {
		testdata->r = testdata->m;
		return;
	}

	// build the cross kernel matrix between train and test
	double *K2 = NULL;
	gensvm_make_crosskernel(model, traindata, testdata, &K2);

	// generate the data matrix N = K2 * M * Sigma^{-2}
	gensvm_make_testfactor(testdata, traindata, K2);

	free(K2);
}

void gensvm_make_kernel(struct GenModel *model, struct GenData *data,
	       	double *K)
{
	long i, j;
	long n = data->n;
	double value;
	double *x1, *x2;

	for (i=0; i<n; i++) {
		for (j=i; j<n; j++) {
			x1 = &data->RAW[i*(data->m+1)+1];
			x2 = &data->RAW[j*(data->m+1)+1];
			if (model->kerneltype == K_POLY)
				value = gensvm_dot_poly(x1, x2,
						model->kernelparam, data->m);
			else if (model->kerneltype == K_RBF)
				value = gensvm_dot_rbf(x1, x2,
						model->kernelparam, data->m);
			else if (model->kerneltype == K_SIGMOID)
				value = gensvm_dot_sigmoid(x1, x2,
						model->kernelparam, data->m);
			else {
				fprintf(stderr, "Unknown kernel type in "
						"gensvm_make_kernel\n");
				exit(1);
			}
			matrix_set(K, n, i, j, value);
			matrix_set(K, n, j, i, value);
		}
	}
}

/**
 * @brief Find the (reduced) eigendecomposition of a kernel matrix.
 *
 * @details.
 * tbd
 *
 */
long gensvm_make_eigen(double *K, long n, double **P, double **Sigma)
{
	int M, status, LWORK, *IWORK, *IFAIL;
	long i, j, num_eigen, cutoff_idx;
	double max_eigen, abstol, *WORK;

	double *tempSigma = Malloc(double, n);
	double *tempP = Malloc(double, n*n);

	IWORK = Malloc(int, 5*n);
	IFAIL = Malloc(int, n);

	// highest precision eigenvalues, may reduce for speed
	abstol = 2.0*dlamch('S');

	// first perform a workspace query to determine optimal size of the
	// WORK array.
	WORK = Malloc(double, 1);
	status = dsyevx(
			'V',
			'A',
			'U',
			n,
			K,
			n,
			0,
			0,
			0,
			0,
			abstol,
			&M,
			tempSigma,
			tempP,
			n,
			WORK,
			-1,
			IWORK,
			IFAIL);
	LWORK = WORK[0];

	// allocate the requested memory for the eigendecomposition
	WORK = (double *)realloc(WORK, LWORK*sizeof(double));
	status = dsyevx(
			'V',
			'A',
			'U',
			n,
			K,
			n,
			0,
			0,
			0,
			0,
			abstol,
			&M,
			tempSigma,
			tempP,
			n,
			WORK,
			LWORK,
			IWORK,
			IFAIL);

	if (status != 0) {
		fprintf(stderr, "Nonzero exit status from dsyevx. Exiting...");
		exit(1);
	}

	// Select the desired number of eigenvalues, depending on their size.
	// dsyevx sorts eigenvalues in ascending order.
	max_eigen = tempSigma[n-1];
	cutoff_idx = 0;

	for (i=0; i<n; i++)
		if (tempSigma[i]/max_eigen > 1e-8 ) {
			cutoff_idx = i;
			break;
		}

	num_eigen = n - cutoff_idx;

	*Sigma = Calloc(double, num_eigen);

	for (i=0; i<num_eigen; i++) {
		(*Sigma)[i] = tempSigma[n-1 - i];
	}

	// revert P to row-major order and copy only the the columns
	// corresponding to the selected eigenvalues
	*P = Calloc(double, n*num_eigen);
	for (j=n-1; j>n-1-num_eigen; j--) {
		for (i=0; i<n; i++) {
			(*P)[i*num_eigen + (n-1)-j] = tempP[i + j*n];
		}
	}

	free(tempSigma);
	free(tempP);
	free(IWORK);
	free(IFAIL);
	free(WORK);

	return num_eigen;
}

void gensvm_make_crosskernel(struct GenModel *model,
	       	struct GenData *data_train, struct GenData *data_test,
	       	double **K2)
{
	long i, j;
	long n_train = data_train->n;
	long n_test = data_test->n;
	long m = data_test->m;
	double value;
	double *x1, *x2;

	*K2 = Calloc(double, n_test*n_train);

	for (i=0; i<n_test; i++) {
		for (j=0; j<n_train; j++) {
			x1 = &data_test->RAW[i*(m+1)+1];
			x2 = &data_train->RAW[j*(m+1)+1];
			if (model->kerneltype == K_POLY)
				value = gensvm_dot_poly(x1, x2,
						model->kernelparam,
					       	m);
			else if (model->kerneltype == K_RBF)
				value = gensvm_dot_rbf(x1, x2,
						model->kernelparam,
					       	m);
			else if (model->kerneltype == K_SIGMOID)
				value = gensvm_dot_sigmoid(x1, x2,
						model->kernelparam,
					       	m);
			else {
				fprintf(stderr, "Unknown kernel type in "
						"gensvm_make_crosskernel\n");
				exit(1);
			}
			matrix_set((*K2), n_train, i, j, value);
		}
	}
}

void gensvm_make_trainfactor(struct GenData *data, double *P, double *Sigma,
	       	long r)
{
	long i, j, n = data->n;
	double value;

	// allocate Z
	data->Z = Calloc(double, n*(r+1));

	// Write data->Z = [1 M] = [1 P*Sigma]
	for (i=0; i<n; i++) {
		for (j=0; j<r; j++) {
			value = matrix_get(P, r, i, j);
			value *= matrix_get(Sigma, 1, j, 0);
			matrix_set(data->Z, r+1, i, j+1, value);
		}
		matrix_set(data->Z, r+1, i, 0, 1.0);
	}

	// Set data->r to r so data knows the width of Z
	data->r = r;
}

void gensvm_make_testfactor(struct GenData *testdata,
	       	struct GenData *traindata, double *K2)
{
	long n1, n2, r, i, j;
	double value,
	       *N = NULL,
	       *M = NULL;

	n1 = traindata->n;
	n2 = testdata->n;
	r = traindata->r;

	N = Calloc(double, n2*r);
	M = Calloc(double, n1*r);

	// copy M from traindata->Z because we need it in dgemm without column
	// of 1's.
	for (i=0; i<n1; i++)
		for (j=0; j<r; j++)
			matrix_set(M, r, i, j,
				       	matrix_get(traindata->Z, r+1, i, j+1));

	// Multiply K2 with M and store in N
	cblas_dgemm(
			CblasRowMajor,
			CblasNoTrans,
			CblasNoTrans,
			n2,
			r,
			n1,
			1.0,
			K2,
			n1,
			M,
			r,
			0.0,
			N,
			r);

	// Multiply N with Sigma^{-2}
	for (j=0; j<r; j++) {
		value = pow(matrix_get(traindata->Sigma, 1, j, 0), -2.0);
		for (i=0; i<n2; i++)
			matrix_mul(N, r, i, j, value);
	}

	// write N to Z with a column of ones
	testdata->Z = Calloc(double, n2*(r+1));
	for (i=0; i<n2; i++) {
		for (j=0; j<r; j++) {
			matrix_set(testdata->Z, r+1, i, j+1,
					matrix_get(N, r, i, j));
		}
		matrix_set(testdata->Z, r+1, i, 0, 1.0);
	}
	// Set r to testdata
	testdata->r = r;

	free(M);
	free(N);
}

/**
 * @brief Compute the RBF kernel between two vectors
 *
 * @details
 * The RBF kernel is computed between two vectors. This kernel is defined as
 * @f[
 * 	k(x_1, x_2) = \exp( -\gamma \| x_1 - x_2 \|^2 )
 * @f]
 * where @f$ \gamma @f$ is a kernel parameter specified.
 *
 * @param[in] 	x1 		first vector
 * @param[in] 	x2 		second vector
 * @param[in] 	kernelparam 	array of kernel parameters (gamma is first
 * 				element)
 * @param[in] 	n 		length of the vectors x1 and x2
 * @returns  			kernel evaluation
 */
double gensvm_dot_rbf(double *x1, double *x2, double *kernelparam, long n)
{
	long i;
	double value = 0.0;

	for (i=0; i<n; i++)
		value += (x1[i] - x2[i]) * (x1[i] - x2[i]);
	value *= -kernelparam[0];
	return exp(value);
}

/**
 * @brief Compute the polynomial kernel between two vectors
 *
 * @details
 * The polynomial kernel is computed between two vectors. This kernel is
 * defined as
 * @f[
 * 	k(x_1, x_2) = ( \gamma \langle x_1, x_2 \rangle + c)^d
 * @f]
 * where @f$ \gamma @f$, @f$ c @f$ and @f$ d @f$ are kernel parameters.
 *
 * @param[in] 	x1 		first vector
 * @param[in] 	x2 		second vector
 * @param[in] 	kernelparam 	array of kernel parameters (gamma, c, d)
 * @param[in] 	n 		length of the vectors x1 and x2
 * @returns 			kernel evaluation
 */
double gensvm_dot_poly(double *x1, double *x2, double *kernelparam, long n)
{
	long i;
	double value = 0.0;
	for (i=0; i<n; i++)
		value += x1[i]*x2[i];
	value *= kernelparam[0];
	value += kernelparam[1];
	return pow(value, ((int) kernelparam[2]));
}

/**
 * @brief Compute the sigmoid kernel between two vectors
 *
 * @details
 * The sigmoid kernel is computed between two vectors. This kernel is defined
 * as
 * @f[
 * 	k(x_1, x_2) = \tanh( \gamma \langle x_1 , x_2 \rangle + c)
 * @f]
 * where @f$ \gamma @f$ and @f$ c @f$ are kernel parameters.
 *
 * @param[in] 	x1 		first vector
 * @param[in] 	x2 		second vector
 * @param[in] 	kernelparam 	array of kernel parameters (gamma, c)
 * @param[in] 	n 		length of the vectors x1 and x2
 * @returns 			kernel evaluation
 */
double gensvm_dot_sigmoid(double *x1, double *x2, double *kernelparam, long n)
{
	long i;
	double value = 0.0;
	for (i=0; i<n; i++)
		value += x1[i]*x2[i];
	value *= kernelparam[0];
	value += kernelparam[1];
	return tanh(value);
}
