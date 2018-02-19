/**
 * @file gensvm_kernel.c
 * @author G.J.J. van den Burg
 * @date 2013-10-18
 * @brief Defines main functions for use of kernels in GenSVM.
 *
 * @details
 * Functions for constructing different kernels using user-supplied
 * parameters. Also contains the functions for decomposing the
 * kernel matrix using several decomposition methods.
 *
 * @copyright
 Copyright 2016, G.J.J. van den Burg.

 This file is part of GenSVM.

 GenSVM is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 GenSVM is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with GenSVM. If not, see <http://www.gnu.org/licenses/>.

*/

#include "gensvm_kernel.h"
#include "gensvm_print.h"

/**
 * @brief Copy the kernelparameters from GenModel to GenData
 *
 * @details
 * This is a little utility function to copy the kernel type and kernel 
 * parameters from a GenModel struct to a GenData struct.
 *
 * @param[in] 	 model 	a GenModel struct
 * @param[in] 	 data 	a GenData struct
 *
 */
void gensvm_kernel_copy_kernelparam_to_data(struct GenModel *model, 
		struct GenData *data)
{
	data->kerneltype = model->kerneltype;
	data->gamma = model->gamma;
	data->coef = model->coef;
	data->degree = model->degree;
}

/**
 * @brief Do the preprocessing steps needed to perform kernel GenSVM
 *
 * @details
 * To achieve nonlinearity through kernels in GenSVM, a preprocessing step is
 * needed. This preprocessing step computes the full kernel matrix, and an
 * eigendecomposition of this matrix. Next, it computes a matrix @f$\textbf{M}
 * = \textbf{P}\boldsymbol{\Sigma}@f$ which takes the role as data matrix in
 * the optimization algorithm.
 *
 * @sa
 * gensvm_kernel_compute(), gensvm_kernel_eigendecomp(), 
 * gensvm_kernel_trainfactor(), gensvm_kernel_postprocess()
 *
 * @param[in] 		model 	input GenSVM model
 * @param[in,out] 	data 	input structure with the data. On exit,
 * 				contains the training factor in GenData::Z,
 * 				and the original data in GenData::RAW
 *
 */
void gensvm_kernel_preprocess(struct GenModel *model, struct GenData *data)
{
	if (model->kerneltype == K_LINEAR) {
		data->r = data->m;
		return;
	}

	long r, n = data->n;
	double *P = NULL,
	       *Sigma = NULL,
	       *K = NULL;

	// build the kernel matrix
	K = Calloc(double, n*n);
	gensvm_kernel_compute(model, data, K);

	// generate the eigen decomposition
	r = gensvm_kernel_eigendecomp(K, n, model->kernel_eigen_cutoff, &P, 
			&Sigma);

	// build M and set to data (leave RAW intact)
	gensvm_kernel_trainfactor(data, P, Sigma, r);

	// Set Sigma to data->Sigma (need it again for prediction)
	if (data->Sigma != NULL) {
		free(data->Sigma);
		data->Sigma = NULL;
	}
	data->Sigma = Sigma;

	// write kernel params to data
	gensvm_kernel_copy_kernelparam_to_data(model, data);

	free(K);
	free(P);
}

/**
 * @brief Compute the kernel postprocessing factor
 *
 * @details
 * This function computes the postprocessing factor needed to do predictions 
 * with kernels in GenSVM. This is a wrapper around gensvm_kernel_cross() and 
 * gensvm_kernel_testfactor().
 *
 * @param[in] 		model 		a GenSVM model
 * @param[in] 		traindata 	the training dataset
 * @param[in,out] 	testdata 	the test dataset. On exit, GenData::Z
 * 					contains the testfactor
 */
void gensvm_kernel_postprocess(struct GenModel *model,
		struct GenData *traindata, struct GenData *testdata)
{
	if (model->kerneltype == K_LINEAR) {
		testdata->r = testdata->m;
		return;
	}

	// build the cross kernel matrix between train and test
	double *K2 = gensvm_kernel_cross(model, traindata, testdata);

	// generate the data matrix N = K2 * M * Sigma^{-2}
	gensvm_kernel_testfactor(testdata, traindata, K2);

	free(K2);
}

/**
 * @brief Compute the kernel matrix
 *
 * @details
 * This function computes the kernel matrix of a data matrix based on the
 * requested kernel type and the kernel parameters. The potential types of
 * kernel functions are document in KernelType. This function uses a naive
 * multiplication and computes the entire upper triangle of the kernel matrix,
 * then copies this over to the lower triangle.
 *
 * @param[in] 	model 	a GenModel structure with the model
 * @param[in] 	data 	a GenData structure with the data
 * @param[out]	K 	an nxn preallocated kernel matrix
 *
 */
void gensvm_kernel_compute(struct GenModel *model, struct GenData *data,
		double *K)
{
	long i, j, incx, incy;
	long n = data->n;
	double value;
	double *x1 = NULL,
	       *x2 = NULL;

	#if MAJOR_ORDER == 'r'
	incx = 1;
	incy = 1;
	#else
	incx = n;
	incy = n;
	#endif

	for (i=0; i<n; i++) {
		for (j=i; j<n; j++) {
			#if MAJOR_ORDER == 'r'
			x1 = &data->RAW[i*(data->m+1)+1];
			x2 = &data->RAW[j*(data->m+1)+1];
			#else
			x1 = &data->RAW[i+n];
			x2 = &data->RAW[j+n];
			#endif
			if (model->kerneltype == K_POLY)
				value = gensvm_kernel_dot_poly(x1, x2, data->m,
						incx, incy, model->gamma, 
						model->coef, model->degree);
			else if (model->kerneltype == K_RBF)
				value = gensvm_kernel_dot_rbf(x1, x2, data->m,
						incx, incy, model-> gamma);
			else if (model->kerneltype == K_SIGMOID)
				value = gensvm_kernel_dot_sigmoid(x1, x2, 
						data->m, incx, incy, 
						model->gamma, model->coef);
			else {
				// LCOV_EXCL_START
				gensvm_error("[GenSVM Error]: Unknown kernel type in "
						"gensvm_make_kernel\n");
				exit(EXIT_FAILURE);
				// LCOV_EXCL_STOP
			}
			matrix_set(K, n, n, i, j, value);
			matrix_set(K, n, n, j, i, value);
		}
	}
}

/**
 * @brief Find the (reduced) eigendecomposition of a kernel matrix
 *
 * @details
 * Compute the reduced eigendecomposition of the kernel matrix. This uses the 
 * LAPACK function dsyevx to do the computation. Only those 
 * eigenvalues/eigenvectors are kept for which the ratio between the 
 * eigenvalue and the largest eigenvalue is larger than cutoff.  This function 
 * uses the highest precision eigenvalues, twice the underflow threshold (see 
 * dsyevx documentation). 
 *
 * @param[in] 		K 		the kernel matrix
 * @param[in] 		n 		the dimension of the kernel matrix
 * @param[in] 		cutoff 		mimimum ratio of eigenvalue to largest
 * 					eigenvalue for the eigenvector to be 
 * 					included
 * @param[out] 		P_ret 		on exit contains the eigenvectors
 * @param[out] 		Sigma_ret 	on exit contains the eigenvalues
 *
 * @return 			the number of eigenvalues kept
 */
long gensvm_kernel_eigendecomp(double *K, long n, double cutoff, double **P_ret,
		double **Sigma_ret)
{
	int M, status, LWORK, *IWORK = NULL,
	    *IFAIL = NULL;
	long i, j, num_eigen, cutoff_idx;
	double max_eigen, abstol, *WORK = NULL,
	       *Sigma = NULL,
	       *P = NULL;

	double *tempSigma = Malloc(double, n);
	double *tempP = Malloc(double, n*n);

	IWORK = Malloc(int, 5*n);
	IFAIL = Malloc(int, n);

	// highest precision eigenvalues, may reduce for speed
	abstol = 2.0*dlamch('S');

	// first perform a workspace query to determine optimal size of the
	// WORK array.
	WORK = Malloc(double, 1);
	status = dsyevx('V', 'A', 'U', n, K, n, 0, 0, 0, 0, abstol, &M,
			tempSigma, tempP, n, WORK, -1, IWORK, IFAIL);
	LWORK = WORK[0];

	// allocate the requested memory for the eigendecomposition
	WORK = (double *)realloc(WORK, LWORK*sizeof(double));
	status = dsyevx('V', 'A', 'U', n, K, n, 0, 0, 0, 0, abstol, &M,
			tempSigma, tempP, n, WORK, LWORK, IWORK, IFAIL);

	if (status != 0) {
		// LCOV_EXCL_START
		gensvm_error("[GenSVM Error]: Nonzero exit status from dsyevx.\n");
		exit(EXIT_FAILURE);
		// LCOV_EXCL_STOP
	}

	// Select the desired number of eigenvalues, depending on their size.
	// dsyevx sorts eigenvalues in ascending order.
	max_eigen = tempSigma[n-1];
	cutoff_idx = 0;

	for (i=0; i<n; i++) {
		if (tempSigma[i]/max_eigen > cutoff) {
			cutoff_idx = i;
			break;
		}
	}

	num_eigen = n - cutoff_idx;

	// In the mathematical derivation (see paper), we state that the 
	// diagonal matrix Sigma contains the square root of the eigenvalues 
	// (i.e. the eigendecomposition is: K = P * Sigma^2 * P').
	Sigma = Calloc(double, num_eigen);
	for (i=0; i<num_eigen; i++) {
		Sigma[i] = sqrt(tempSigma[n-1 - i]);
	}

	// revert P to row-major order and copy only the the columns
	// corresponding to the selected eigenvalues
	P = Calloc(double, n*num_eigen);
	#if MAJOR_ORDER == 'r'
	for (j=n-1; j>n-1-num_eigen; j--) {
		for (i=0; i<n; i++) {
			P[i*num_eigen + (n-1)-j] = tempP[i + j*n];
		}
	}
	#else
	for (j=n-1; j>n-1-num_eigen; j--) {
		for (i=0; i<n; i++) {
			P[i + (n - 1 - j)*n] = tempP[i + j*n];
		}
	}
	#endif

	free(tempSigma);
	free(tempP);
	free(IWORK);
	free(IFAIL);
	free(WORK);

	*Sigma_ret = Sigma;
	*P_ret = P;

	return num_eigen;
}

/**
 * @brief Compute the kernel crossproduct between two datasets
 *
 * @details
 * Given a training set @f$\textbf{X}@f$ with feature space mapping 
 * @f$\boldsymbol{\Phi}@f$ and a test set @f$\textbf{X}_2@f$ with feature 
 * space mapping @f$\boldsymbol{\Phi}_2@f$, the crosskernel @f$\textbf{K}_2@f$ 
 * is given by @f$\textbf{K}_2 = \boldsymbol{\Phi}_2 \boldsymbol{\Phi}'@f$.  
 * Thus, an element in row @f$i@f$ and column @f$j@f$ in @f$\textbf{K}_2@f$ 
 * equals the kernel product between the @f$i@f$-th row of @f$\textbf{X}_2@f$ 
 * and the @f$j@f$-th row of @f$\textbf{X}@f$.
 *
 * @param[in] 	model 		the GenSVM model
 * @param[in] 	data_train 	the training dataset
 * @param[in] 	data_test 	the test dataset
 *
 * @return  	the matrix @f$\textbf{K}_2@f$
 */
double *gensvm_kernel_cross(struct GenModel *model, struct GenData *data_train,
		struct GenData *data_test)
{
	long i, j, incx, incy;
	long n_train = data_train->n;
	long n_test = data_test->n;
	long m = data_test->m;
	double value, *x1 = NULL,
	       *x2 = NULL,
	       *K2 = Calloc(double, n_test * n_train);

	#if MAJOR_ORDER == 'r'
	incx = 1;
	incy = 1;
	#else
	incx = n_test;
	incy = n_train;
	#endif

	for (i=0; i<n_test; i++) {
		for (j=0; j<n_train; j++) {
			#if MAJOR_ORDER == 'r'
			x1 = &data_test->RAW[i*(m+1)+1];
			x2 = &data_train->RAW[j*(m+1)+1];
			#else
			x1 = &data_test->RAW[i+n_test];
			x2 = &data_train->RAW[j+n_train];
			#endif
			if (model->kerneltype == K_POLY)
				value = gensvm_kernel_dot_poly(x1, x2, m, 
						incx, incy, model->gamma, 
						model->coef, model->degree);
			else if (model->kerneltype == K_RBF)
				value = gensvm_kernel_dot_rbf(x1, x2, m,
						incx, incy, model->gamma);
			else if (model->kerneltype == K_SIGMOID)
				value = gensvm_kernel_dot_sigmoid(x1, x2, m,
						incx, incy, model->gamma, 
						model->coef);
			else {
				// LCOV_EXCL_START
				gensvm_error("[GenSVM Error]: Unknown kernel type in "
						"gensvm_make_crosskernel\n");
				exit(EXIT_FAILURE);
				// LCOV_EXCL_STOP
			}
			matrix_set(K2, n_test, n_train, i, j, value);
		}
	}
	return K2;
}

/**
 * @brief Compute the training factor as part of kernel preprocessing
 *
 * @details
 * This function computes the matrix product @f$\textbf{M} = 
 * \textbf{P}\boldsymbol{\Sigma}@f$ and stores the result in GenData::Z, 
 * preceded by a column of ones. It also sets GenData::r to the number of 
 * eigenvectors that were includedin P and Sigma. Note that P and Sigma 
 * correspond to the reduced eigendecomposition of the kernel matrix.
 *
 * @param[in,out] 	data 	a GenData structure. On exit, GenData::Z and
 * 				GenData::r are updated as described above.
 * @param[in] 		P 	the eigenvectors
 * @param[in] 		Sigma 	the eigenvalues
 * @param[in] 		r 	the number of eigenvalues and eigenvectors
 */
void gensvm_kernel_trainfactor(struct GenData *data, double *P, double *Sigma,
		long r)
{
	long i, j, n = data->n;
	double value;

	// allocate Z
	data->Z = Calloc(double, n*(r+1));

	// Write data->Z = [1 M] = [1 P*Sigma]
	for (i=0; i<n; i++) {
		for (j=0; j<r; j++) {
			value = matrix_get(P, n, r, i, j);
			value *= matrix_get(Sigma, r, 1, j, 0);
			matrix_set(data->Z, n, r+1, i, j+1, value);
		}
		matrix_set(data->Z, n, r+1, i, 0, 1.0);
	}

	// Set data->r to r so data knows the width of Z
	data->r = r;
}

/**
 * @brief Calculate the matrix product for the testfactor
 *
 * @details
 * To predict class labels when kernels are used, a transformation of the
 * testdata has to be performed to get the simplex space vectors. This
 * transformation is based on the matrix @f$\textbf{K}_2@f$ (as calculated by
 * gensvm_make_crosskernel()) and the matrices @f$\textbf{M} = 
 * \textbf{P}*\boldsymbol{\Sigma}@f$) and @f$\boldsymbol{\Sigma}@f$. The 
 * testfactor is equal to @f$\textbf{K}_2 \textbf{M} 
 * \boldsymbol{\Sigma}^{-2}@f$.
 *
 * @param[out] 	testdata 	a GenData struct with the testdata, contains
 * 				the testfactor in GenData::Z on exit preceded 
 * 				by a column of ones.
 * @param[in] 	traindata 	a GenData struct with the training data
 * @param[in] 	K2 		crosskernel between the train and test data
 */
void gensvm_kernel_testfactor(struct GenData *testdata,
		struct GenData *traindata, double *K2)
{
	long n1, n2, r, i, j;
	double value, *N = NULL,
	       *M = NULL;

	n1 = traindata->n;
	n2 = testdata->n;
	r = traindata->r;

	N = Calloc(double, n2*r);
	M = Calloc(double, n1*r);

	// copy M from traindata->Z because we need it in dgemm without column
	// of 1's.
	for (i=0; i<n1; i++) {
		for (j=0; j<r; j++) {
			value = matrix_get(traindata->Z, n1, r+1, i, j+1);
			matrix_set(M, n1, r, i, j, value);
		}
	}

	// Multiply K2 with M and store in N
	#if MAJOR_ORDER == 'r'
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n2, r, n1, 1.0,
			K2, n1, M, r, 0.0, N, r);
	#else
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n2, r, n1, 1.0,
			K2, n2, M, n1, 0.0, N, n2);
	#endif

	// Multiply N with Sigma^{-2}
	for (j=0; j<r; j++) {
		value = pow(matrix_get(traindata->Sigma, r, 1, j, 0), -2.0);
		for (i=0; i<n2; i++)
			matrix_mul(N, n2, r, i, j, value);
	}

	// write N to Z with a column of ones
	testdata->Z = Calloc(double, n2*(r+1));
	for (i=0; i<n2; i++) {
		for (j=0; j<r; j++) {
			matrix_set(testdata->Z, n2, r+1, i, j+1,
					matrix_get(N, n2, r, i, j));
		}
		matrix_set(testdata->Z, n2, r+1, i, 0, 1.0);
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
 * @param[in] 	n 		length of the vectors x1 and x2
 * @param[in] 	gamma 		gamma parameter of the kernel
 * @returns  			kernel evaluation
 */
double gensvm_kernel_dot_rbf(double *x, double *y, long n, 
		long incx, long incy, double gamma)
{
	long i, ix, iy;
	double value = 0.0;

	if (incx == 1 && incy == 1) {
		for (i=0; i<n; i++) {
			value += (x[i] - y[i]) * (x[i] - y[i]);
		}
	} else {
		ix = 0;
		iy = 0;
		for (i=0; i<n; i++) {
			value += (x[ix] - y[iy]) * (x[ix] - y[iy]);
			ix += incx;
			iy += incy;
		}
	}

	value *= -gamma;
	return exp(value);
}

/**
 * @brief Compute the polynomial kernel between two vectors
 *
 * @details
 * The polynomial kernel is computed between two vectors. This kernel is
 * defined as
 * @f[
 * 	k(x_1, x_2) = ( \gamma \langle x_1, x_2 \rangle + coef)^{degree}
 * @f]
 * where @f$ \gamma @f$, @f$ coef @f$ and @f$ degree @f$ are kernel 
 * parameters.
 *
 * @param[in] 	x1 		first vector
 * @param[in] 	x2 		second vector
 * @param[in] 	n 		length of the vectors x1 and x2
 * @param[in] 	gamma 		gamma parameter of the kernel
 * @param[in] 	coef 		coef parameter of the kernel
 * @param[in] 	degree 		degree parameter of the kernel
 * @returns 			kernel evaluation
 */
double gensvm_kernel_dot_poly(double *x, double *y, long n, long incx, 
		long incy, double gamma, double coef, double degree)
{
	double value = cblas_ddot(n, x, incx, y, incy);
	value *= gamma;
	value += coef;
	return pow(value, degree);
}

/**
 * @brief Compute the sigmoid kernel between two vectors
 *
 * @details
 * The sigmoid kernel is computed between two vectors. This kernel is defined
 * as
 * @f[
 * 	k(x_1, x_2) = \tanh( \gamma \langle x_1 , x_2 \rangle + coef)
 * @f]
 * where @f$ \gamma @f$ and @f$ coef @f$ are kernel parameters.
 *
 * @param[in] 	x1 		first vector
 * @param[in] 	x2 		second vector
 * @param[in] 	n 		length of the vectors x1 and x2
 * @param[in] 	gamma 		gamma parameter of the kernel
 * @param[in] 	coef 		coef parameter of the kernel
 * @returns 			kernel evaluation
 */
double gensvm_kernel_dot_sigmoid(double *x, double *y, long n, long incx, 
		long incy, double gamma, double coef)
{
	double value = cblas_ddot(n, x, incx, y, incy);
	value *= gamma;
	value += coef;
	return tanh(value);
}

/**
 * @brief Compute the eigenvalues and optionally the eigenvectors of a
 * symmetric matrix.
 *
 * @details
 * This is a wrapper function around the external LAPACK function.
 *
 * See the LAPACK documentation at:
 * http://www.netlib.org/lapack/explore-html/d2/d97/dsyevx_8f.html
 *
 */
int dsyevx(char JOBZ, char RANGE, char UPLO, int N, double *A, int LDA,
		double VL, double VU, int IL, int IU, double ABSTOL, int *M,
		double *W, double *Z, int LDZ, double *WORK, int LWORK,
		int *IWORK, int *IFAIL)
{
	extern void dsyevx_(char *JOBZ, char *RANGE, char *UPLO, int *Np,
			double *A, int *LDAp, double *VLp, double *VUp,
			int *ILp, int *IUp, double *ABSTOLp, int *M,
			double *W, double *Z, int *LDZp, double *WORK,
			int *LWORKp, int *IWORK, int *IFAIL, int *INFOp);
	int INFO;
	dsyevx_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &IU, &ABSTOL,
			M, W, Z, &LDZ, WORK, &LWORK, IWORK, IFAIL, &INFO);
	return INFO;
}

/**
 * @brief Determine double precision machine parameters.
 *
 * @details
 * This is a wrapper function around the external LAPACK function.
 *
 * See the LAPACK documentation at:
 * http://www.netlib.org/lapack/explore-html/d5/dd4/dlamch_8f.html
 */
double dlamch(char CMACH)
{
	extern double dlamch_(char *CMACH);
	return dlamch_(&CMACH);
}
