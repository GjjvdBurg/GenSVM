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
	if (model->kerneltype == K_LINEAR)
		return;

	long i, j;
	long n = model->n;
	double value;
	double *x1, *x2;
	double *K = Calloc(double, n*n*sizeof(double));

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
				value = msvmmaj_compute_rbf(x1, x2, 
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

	print_matrix(K, n, n);

	double *P = Malloc(double, n*n);
	double *Lambda = Malloc(double, n);
	long num_eigen = msvmmaj_make_eigen(K, n, P, Lambda);

	// copy eigendecomp to data
	data->Z = realloc(data->Z, n*(n+1)*sizeof(double));
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++)
			matrix_set(data->Z, n+1, i, j+1,
					matrix_get(P, n, i, j));
		matrix_set(data->Z, n+1, i, 0, 1.0);
	}
	data->m = n;

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
	model->m = n;
	free(K);
}

/**
 * @brief Find the eigendecomposition of a kernel matrix.
 *
 * @details.
 * tbd
 *
 *
 */
long msvmmaj_make_eigen(double *K, long n, double *P, double *Lambda)
{
	int M, status, LWORK, *IWORK, *IFAIL;
	double abstol, *WORK;

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
			Lambda,
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
			Lambda,
			tempP,
			n,
			WORK,
			LWORK,
			IWORK,
			IFAIL);

	printf("status = %i\n", status);
	printf("Number of eigenvalues found: %i\n", M);

	if (status != 0) {
		fprintf(stderr, "Nonzero exit status from dsyevx. Exiting...");
		exit(1);
	}
	// Here you can put code to select the necessary eigenvectors, 
	// depending on the size of the eigenvalues.
	// For now, let's just print the eigenvalues and exit
	
	print_matrix(Lambda, n, 1);
	
	// revert P to row-major order	
	long i, j;
	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
			P[i*n+j] = tempP[j*n+i];

	print_matrix(P, n, n);

	free(tempP);

	// replace by number of columns of P
	return n;
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
