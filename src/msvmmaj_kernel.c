/**
 * @file msvmmaj_kernel.c
 * @author Gertjan van den Burg
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

#include "msvmmaj.h"
#include "msvmmaj_kernel.h"
#include "msvmmaj_lapack.h"
#include "msvmmaj_matrix.h"
#include "util.h"

/**
 * @brief Create the kernel matrix
 *
 * Create a kernel matrix based on the specified kerneltype. Kernel parameters 
 * are assumed to be specified in the model.
 *
 * @param[in] 	model 	MajModel specifying the parameters
 * @param[in] 	data 	MajData specifying the data.
 *
 */
void msvmmaj_make_kernel(struct MajModel *model, struct MajData *data)
{
	long i, j;
	// Determine if a kernel needs to be computed. This is not the case if 
	// a LINEAR kernel is requested in the model, or if the requested 
	// kernel is already in the data. 

	if (model->kerneltype == K_LINEAR) {
		data->J = Calloc(double, data->m+1);
		for (i=1; i<data->m+1; i++) {
			matrix_set(data->J, 1, i, 0, 1.0);
		}
		return;
	}

	/*
	switch (model->kerneltype) {
		case K_LINEAR:
			// if data has another kernel, free that matrix and 
			// assign Z to RAW
			if (data->kerneltype != K_LINEAR) {
				free(data->Z);
				data->Z = data->RAW;
			}
			data->J = Calloc(double, data->m+1);
			for (i=1; i<model->m+1; i++) {
				matrix_set(data->J, 1, i, 0, 1.0);
			}
			return;
		case K_POLY:
			// if data has another kernel, we need to recalculate
			if (data->kerneltype != K_POLY) {
				break;
			}
			// if it is poly, we only recalculate if the kernel 
			// parameters differ
			if (data->kernelparam[0] == model->kernelparam[0] &&
				data->kernelparam[1] == model->kernelparam[1] &&
				data->kernelparam[2] == model->kernelparam[2])
				// < do something with J ?
				return;
		case K_RBF:
			if (data->kerneltype != K_RBF)
				break;
			if (data->kernelparam[0] == model->kernelparam[0])
				// < do something with J ?
				return;
		case K_SIGMOID:
			if (data->kerneltype != K_SIGMOID)
				break;
			if (data->kernelparam[0] == model->kernelparam[0] &&
				data->kernelparam[1] == model->kernelparam[1])
				// < do something with J ?
				return;
	}
	*/
	long n = data->n;
	double value;
	double *x1, *x2;
	double *K = Calloc(double, n*n);

	for (i=0; i<n; i++) {
		for (j=i; j<n; j++) {
			x1 = &data->Z[i*(data->m+1)+1];
			x2 = &data->Z[j*(data->m+1)+1];
			if (model->kerneltype == K_POLY)
				value = msvmmaj_compute_poly(x1, x2, 
						model->kernelparam, data->m);
			else if (model->kerneltype == K_RBF)
				value = msvmmaj_compute_rbf(x1, x2, 
						model->kernelparam, data->m);
			else if (model->kerneltype == K_SIGMOID)
				value = msvmmaj_compute_sigmoid(x1, x2, 
						model->kernelparam, data->m);
			else {
				fprintf(stderr, "Unknown kernel type in "
						"msvmmaj_make_kernel\n");
				exit(1);
			}
			matrix_set(K, n, i, j, value);
			matrix_set(K, n, j, i, value);
		}
	}

	double *P = NULL;
	double *Sigma = NULL;
	long num_eigen = msvmmaj_make_eigen(K, n, &P, &Sigma);
	//printf("num eigen: %li\n", num_eigen);
	data->m = num_eigen;

	// copy eigendecomp to data
	data->Z = Calloc(double, n*(num_eigen+1));
	for (i=0; i<n; i++) {
		for (j=0; j<num_eigen; j++) {
			value = matrix_get(P, num_eigen, i, j);
			matrix_set(data->Z, num_eigen+1, i, j, value);
		}
		matrix_set(data->Z, num_eigen+1, i, 0, 1.0);
	}

	// Set the regularization matrix (change if not full rank used)
	data->J = Calloc(double, data->m+1);
	for (i=1; i<data->m+1; i++) {
		value = 1.0/matrix_get(Sigma, 1, i-1, 0);
		matrix_set(data->J, 1, i, 0, value);
	}

	// let data know what it's made of
	data->kerneltype = model->kerneltype;
	free(data->kernelparam);
	switch (model->kerneltype) {
		case K_LINEAR:
			break;
		case K_POLY:
			data->kernelparam = Calloc(double, 3);
			data->kernelparam[0] = model->kernelparam[0];
			data->kernelparam[1] = model->kernelparam[1];
			data->kernelparam[2] = model->kernelparam[2];
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
	free(Sigma);
	free(P);
}

/**
 * @brief Find the (reduced) eigendecomposition of a kernel matrix.
 *
 * @details.
 * tbd
 *
 */
long msvmmaj_make_eigen(double *K, long n, double **P, double **Sigma)
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
	//
	max_eigen = tempSigma[n-1];
	cutoff_idx = 0;

	for (i=0; i<n; i++)
		if (tempSigma[i]/max_eigen > 1e-10 ) {
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
	//
	*P = Calloc(double, n*num_eigen);	
	for (j=n-1; j>n-1-num_eigen; j--) {
		for (i=0; i<n; i++) {
			(*P)[i*num_eigen + (n-1)-j] = tempP[i + j*n];
		}
	}

	free(tempSigma);	
	free(tempP);

	return num_eigen;
}

void msvmmaj_make_crosskernel(struct MajModel *model,
	       	struct MajData *data_train, struct MajData *data_test,
	       	double **K2)
{
	long i, j;
	long n_train = data_train->n;
	long n_test = data_test->n;
	long m = data_test->m;
	double value;
	double *x1, *x2;

	*K2 = Calloc(double, n_test*n_train);

	//printf("Training RAW\n");
	//print_matrix(data_train->RAW, n_train, m+1);

	//printf("Testing RAW\n");
	//print_matrix(data_test->RAW, n_test, m+1);

	for (i=0; i<n_test; i++) {
		for (j=0; j<n_train; j++) {
			x1 = &data_test->RAW[i*(m+1)+1];
			x2 = &data_train->RAW[j*(m+1)+1];
			if (model->kerneltype == K_POLY)
				value = msvmmaj_compute_poly(x1, x2,
						model->kernelparam,
					       	m);
			else if (model->kerneltype == K_RBF)
				value = msvmmaj_compute_rbf(x1, x2,
						model->kernelparam,
					       	m);
			else if (model->kerneltype == K_SIGMOID)
				value = msvmmaj_compute_sigmoid(x1, x2,
						model->kernelparam,
					       	m);
			else {
				fprintf(stderr, "Unknown kernel type in "
						"msvmmaj_make_crosskernel\n");
				exit(1);
			}
			matrix_set((*K2), n_train, i, j, value);
		}
	}

	//printf("cross K2:\n");
	//print_matrix((*K2), n_test, n_train);

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
double msvmmaj_compute_poly(double *x1, double *x2, double *kernelparam, long n)
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
double msvmmaj_compute_sigmoid(double *x1, double *x2, double *kernelparam, long n)
{
	long i;
	double value = 0.0;
	for (i=0; i<n; i++)
		value += x1[i]*x2[i];
	value *= kernelparam[0];
	value += kernelparam[1];
	return tanh(value);
}
