/**
 * @file msvmmaj_train.c
 * @author Gertjan van den Burg
 * @date August 9, 2013
 * @brief Main functions for training the MSVMMaj solution.
 * 
 * @details
 * Contains update and loss functions used to actually find
 * the optimal V.
 *
 */

#include <math.h>
#include <cblas.h>

#include "libMSVMMaj.h"
#include "msvmmaj.h"
#include "msvmmaj_lapack.h"
#include "msvmmaj_matrix.h"
#include "msvmmaj_sv.h"
#include "msvmmaj_train.h"
#include "util.h"

/**
 * Maximum number of iterations of the algorithm.
 */
#define MAX_ITER 1000000000

/**
 * @brief The main training loop for MSVMMaj
 *
 * @details
 * This function is the main training function. This function
 * handles the optimization of the model with the given model parameters, with
 * the data given. On return the matrix MajModel::V contains the optimal 
 * weight matrix.
 *
 * In this function, step doubling is used in the majorization algorithm after
 * a burn-in of 50 iterations. If the training is finished, MajModel::t and 
 * MajModel::W are extracted from MajModel::V.
 *
 * @param[in,out] 	model 	the MajModel to be trained. Contains optimal 
 * 				V on exit.
 * @param[in] 		data 	the MajData to train the model with.
 */
void msvmmaj_optimize(struct MajModel *model, struct MajData *data)
{
	long i, j, it = 0;
	double L, Lbar;

	long n = model->n;
	long m = model->m;
	long K = model->K;

	double *B = Calloc(double, n*(K-1));
	double *ZV = Calloc(double, n*(K-1));
	double *ZAZ = Calloc(double, (m+1)*(m+1));	
	double *ZAZV = Calloc(double, (m+1)*(K-1));
	double *ZAZVT = Calloc(double, (m+1)*(K-1));

	note("Starting main loop.\n");
	note("Dataset:\n");
	note("\tn = %i\n", n);
	note("\tm = %i\n", m);
	note("\tK = %i\n", K);
	note("Parameters:\n");
	note("\tkappa = %f\n", model->kappa);
	note("\tp = %f\n", model->p);
	note("\tlambda = %15.16f\n", model->lambda);
	note("\tepsilon = %g\n", model->epsilon);
	note("\n");

	msvmmaj_simplex_gen(model->K, model->U);
	msvmmaj_simplex_diff(model, data);
	msvmmaj_category_matrix(model, data);

	L = msvmmaj_get_loss(model, data, ZV);
	Lbar = L + 2.0*model->epsilon*L;

	printf("V:\n");
	print_matrix(model->V, m+1, K-1);

	printf("Q:\n");
	print_matrix(model->Q, n, K);

	printf("H:\n");
	print_matrix(model->H, n, K);

	while ((it < MAX_ITER) && (Lbar - L)/L > model->epsilon)
	{
		// ensure V contains newest V and Vbar contains V from 
		// previous
		msvmmaj_get_update(model, data, B, ZAZ, ZAZV, ZAZVT);
		if (it > 50)
			msvmmaj_step_doubling(model);
		
		Lbar = L;
		L = msvmmaj_get_loss(model, data, ZV);

		if (it%1 == 0) 
			note("iter = %li, L = %15.16f, Lbar = %15.16f, "
			     "reldiff = %15.16f\n", it, L, Lbar, (Lbar - L)/L);
		it++;
	}
	if (L > Lbar)
		fprintf(stderr, "Negative step occurred in majorization.\n");

	note("optimization finished, iter = %li, loss = %15.16f, "
			"rel. diff. = %15.16f\n", it-1, L,
			(Lbar - L)/L);
	note("number of support vectors: %li\n", msvmmaj_num_sv(model, data));

	model->training_error = (Lbar - L)/L;

	for (i=0; i<K-1; i++)
		model->t[i] = matrix_get(model->V, K-1, 0, i);
	for (i=1; i<m+1; i++)
		for (j=0; j<K-1; j++)
			matrix_set(model->W, K-1, i-1, j, 
					matrix_get(model->V, K-1, i, j));
	free(B);
	free(ZV);
	free(ZAZ);
	free(ZAZV);
	free(ZAZVT);
}

/**
 * @brief Calculate the current value of the loss function
 * 
 * @details
 * The current loss function value is calculated based on the matrix V in the 
 * given model. Note that the matrix ZV is passed explicitly to avoid having
 * to reallocate memory at every step.
 *
 * @param[in] 		model 	MajModel structure which holds the current 
 * 				estimate V
 * @param[in]  		data 	MajData structure
 * @param[in,out] 	ZV 	pre-allocated matrix ZV which is updated on 
 * 				output
 * @returns 			the current value of the loss function
 */
double msvmmaj_get_loss(struct MajModel *model, struct MajData *data, 
		double *ZV)
{
	long i, j;
	long n = data->n;
	long K = data->K;
	long m = data->m;

	double value, rowvalue, loss = 0.0;

	msvmmaj_calculate_errors(model, data, ZV);
	msvmmaj_calculate_huber(model);

	for (i=0; i<n; i++) {
		rowvalue = 0;
		value = 0;
		for (j=0; j<K; j++) {
			value = matrix_get(model->H, K, i, j);
			value = pow(value, model->p);
			value *= matrix_get(model->R, K, i, j);
			rowvalue += value;
		}
		rowvalue = pow(rowvalue, 1.0/(model->p));
		rowvalue *= model->rho[i];
		loss += rowvalue;
	}
	loss /= ((double) n);

	value = 0;
	for (i=0; i<m+1; i++) {
		rowvalue = 0;
		for (j=0; j<K-1; j++) {
			rowvalue += pow(matrix_get(model->V, K-1, i, j), 2.0);
		}
		value += data->J[i] * rowvalue;
	}
	loss += model->lambda * value;

	return loss;
}

/**
 * @brief Perform a single step of the majorization algorithm to update V
 *
 * @details
 * This function contains the main update calculations of the algorithm. These
 * calculations are necessary to find a new update V. The calculations exist of
 * recalculating the majorization coefficients for all instances and all
 * classes, and solving a linear system to find V.
 *
 * Because the function msvmmaj_get_update() is always called after a call to
 * msvmmaj_get_loss() with the same MajModel::V, it is unnecessary to calculate
 * the updated errors MajModel::Q and MajModel::H here too. This saves on 
 * computation time.
 *
 * In calculating the majorization coefficients we calculate the elements of a 
 * diagonal matrix A with elements
 * @f[
 * 	A_{i, i} = \frac{1}{n} \rho_i \sum_{j \neq k} \left[ 
 * 		\varepsilon_i a_{ijk}^{(p)} + (1 - \varepsilon_i) \omega_i
 * 		a_{ijk}^{(p)} \right],
 * @f]
 * where @f$ k = y_i @f$.
 * Since this matrix is only used to calculate the matrix @f$ Z' A Z @f$, it is
 * efficient to update a matrix ZAZ through consecutive rank 1 updates with 
 * a single element of A and the corresponding row of Z. The BLAS function
 * dsyr is used for this.
 *
 * The B matrix is has rows
 * @f[
 * 	\boldsymbol{\beta}_i' = \frac{1}{n} \rho_i \sum_{j \neq k} \left[
 * 		\varepsilon_i \left( b_{ijk}^{(1)} - a_{ijk}^{(1)}
 * 			\overline{q}_i^{(kj)} \right) + (1 - \varepsilon_i)
 * 		\omega_i \left( b_{ijk}^{(p)} - a_{ijk}^{(p)}
 * 			\overline{q}_i^{(kj)} \right) \right]
 * 		\boldsymbol{\delta}_{kj}'
 * @f]
 * This is also split into two cases, one for which @f$ \varepsilon_i = 1 @f$,
 * and one for when it is 0. The 3D simplex difference matrix is used here, in
 * the form of the @f$ \boldsymbol{\delta}_{kj}' @f$.
 *
 * Finally, the following system is solved
 * @f[
 * 	(\textbf{Z}'\textbf{AZ} + \lambda \textbf{J})\textbf{V} =
 * 		(\textbf{Z}'\textbf{AZ}\overline{\textbf{V}} + \textbf{Z}'
 * 		\textbf{B})
 * @f]
 * solving this system is done through dposv().
 *
 * @todo
 * Consider allocating IPIV and WORK at a higher level, they probably don't 
 * change much during the iterations.
 *
 * @param [in,out] 	model 	model to be updated
 * @param [in] 		data 	data used in model
 * @param [in] 		B 	pre-allocated matrix used for linear coefficients
 * @param [in] 		ZAZ 	pre-allocated matrix used in system
 * @param [in] 		ZAZV 	pre-allocated matrix used in system solving
 * @param [in] 		ZAZVT 	pre-allocated matrix used in system solving
 */
void msvmmaj_get_update(struct MajModel *model, struct MajData *data, double *B,
		double *ZAZ, double *ZAZV, double *ZAZVT)
{
	int status, class;
	long i, j, k;
	double Avalue, Bvalue;
	double omega, value, a, b, q, h, r;

	long n = model->n;
	long m = model->m;
	long K = model->K;

	double kappa = model->kappa;
	double p = model->p;
	double *rho = model->rho;

	// constants which are used often throughout
	const double a2g2 = 0.25*p*(2.0*p - 1.0)*pow((kappa+1.0)/2.0,p-2.0);
	const double in = 1.0/((double) n);

	// clear matrices
	Memset(B, double, n*(K-1));
	Memset(ZAZ, double, (m+1)*(m+1));

	b = 0;
	for (i=0; i<n; i++) {
		value = 0;
		omega = 0;
		for (j=0; j<K; j++) {
			h = matrix_get(model->H, K, i, j);
			r = matrix_get(model->R, K, i, j);
			value += (h*r > 0) ? 1 : 0;
			omega += pow(h, p)*r;
		}
		class = (value <= 1.0) ? 1 : 0;
		omega = (1.0/p)*pow(omega, 1.0/p - 1.0);

		Avalue = 0;
		if (class == 1) {
			for (j=0; j<K; j++) {
				q = matrix_get(model->Q, K, i, j);
				if (q <= -kappa) {
					a = 0.25/(0.5 - kappa/2.0 - q);
					b = 0.5;
				} else if (q <= 1.0) {
					a = 1.0/(2.0*kappa + 2.0);
					b = (1.0 - q)*a;
				} else {
					a = -0.25/(0.5 - kappa/2.0 - q);
					b = 0;
				}
				for (k=0; k<K-1; k++) {
					Bvalue = in*rho[i]*b*matrix3_get(
						model->UU, K-1, K, i, k, j);
					matrix_add(B, K-1, i, k, Bvalue);
				}
				Avalue += a*matrix_get(model->R, K, i, j);
			}
		} else {
			if (2.0 - p < 0.0001) {
				for (j=0; j<K; j++) {
					q = matrix_get(model->Q, K, i, j);
					if (q <= -kappa) {
						b = 0.5 - kappa/2.0 - q;
					} else if ( q <= 1.0) {
						b = pow(1.0 - q, 3.0)/(
							2.0*pow(kappa + 1.0, 
								2.0));
					} else {
						b = 0;
					}
					for (k=0; k<K-1; k++) {
						Bvalue = in*rho[i]*omega*b*
							matrix3_get(
								model->UU,
								K-1,
							       	K,
								i,
								k,
								j);
						matrix_add(
								B,
								K-1,
								i,
								k,
								Bvalue);
					}
				}
				Avalue = 1.5*(K - 1.0);
			} else {
				for (j=0; j<K; j++) {
					q = matrix_get(model->Q, K, i, j);
					if (q <= (p + kappa - 1.0)/(p - 2.0)) {
						a = 0.25*pow(p, 2.0)*pow(
							0.5 - kappa/2.0 - q,
						       		p - 2.0);
					} else if (q <= 1.0) {
						a = a2g2;
					} else {
						a = 0.25*pow(p, 2.0)*pow(
							(p/(p - 2.0))*
							(0.5 - kappa/2.0 - q), 
							p - 2.0);
						b = a*(2.0*q + kappa - 1.0)/
							(p - 2.0) + 
							0.5*p*pow(
								p/(p - 2.0)*
								(0.5 - kappa/
								 2.0 - q), 
								p - 1.0);
					}
					if (q <= -kappa) {
						b = 0.5*p*pow(
							0.5 - kappa/2.0 - q, 
							p - 1.0);
					} else if ( q <= 1.0) {
						b = p*pow(1.0 - q, 
								2.0*p - 1.0)/
							pow(2*kappa+2.0, p);
					}
					for (k=0; k<K-1; k++) {
						Bvalue = in*rho[i]*omega*b*
							matrix3_get(
								model->UU,
								K-1,
								K,
								i,
								k,
								j);
						matrix_add(
								B,
								K-1,
								i,
								k,
								Bvalue);
					}
					Avalue += a*matrix_get(model->R,
						       	K, i, j);
				}
			}
			Avalue *= omega;
		}
		Avalue *= in * rho[i];

		// Now we calculate the matrix ZAZ. Since this is 
		// guaranteed to be symmetric, we only calculate the 
		// upper part of the matrix, and then copy this over
		// to the lower part after all calculations are done.
		// Note that the use of dsym is faster than dspr, even
		// though dspr uses less memory.
		cblas_dsyr(
				CblasRowMajor,
				CblasUpper,
				m+1,
				Avalue,
				&data->Z[i*(m+1)],
				1,
				ZAZ,
				m+1);
	}
	// Copy upper to lower (necessary because we need to switch
	// to Col-Major order for LAPACK).
	/*
	for (i=0; i<m+1; i++)
		for (j=0; j<m+1; j++)
			matrix_set(ZAZ, m+1, j, i, matrix_get(ZAZ, m+1, i, j));
	*/
	
	// Calculate the right hand side of the system we 
	// want to solve.
	cblas_dsymm(
			CblasRowMajor,
			CblasLeft,
			CblasUpper,
			m+1,
			K-1,
			1.0,
			ZAZ,
			m+1,
			model->V,
			K-1,
			0.0,
			ZAZV, 
			K-1);

	cblas_dgemm(
			CblasRowMajor,
			CblasTrans,
			CblasNoTrans,
			m+1,
			K-1,
			n,
			1.0,
			data->Z,
			m+1,
			B,
			K-1,
			1.0,
			ZAZV,
			K-1);

	/* 
	 * Add lambda to all diagonal elements except the first one. Recall 
	 * that ZAZ is of size m+1 and is symmetric.
	 */
	i = 0;
	for (j=0; j<m; j++) {
		i += (m+1) + 1;
		ZAZ[i] += model->lambda * data->J[j+1];
	}
	
	// For the LAPACK call we need to switch to Column-
	// Major order. This is unnecessary for the matrix
	// ZAZ because it is symmetric. The matrix ZAZV 
	// must be converted however.
	for (i=0; i<m+1; i++)
		for (j=0; j<K-1; j++)
			ZAZVT[j*(m+1)+i] = ZAZV[i*(K-1)+j];
	
	// We use the lower ('L') part of the matrix ZAZ, 
	// because we have used the upper part in the BLAS
	// calls above in Row-major order, and Lapack uses
	// column major order. 

	status = dposv(
			'L',
			m+1,
			K-1,
			ZAZ,
			m+1,
			ZAZVT,
			m+1);

	if (status != 0) {
		// This step should not be necessary, as the matrix
		// ZAZ is positive semi-definite by definition. It 
		// is included for safety.
		fprintf(stderr, "Received nonzero status from dposv: %i\n", 
				status);
		int *IPIV = malloc((m+1)*sizeof(int));
		double *WORK = malloc(1*sizeof(double));
		status = dsysv(
				'L',
				m+1,
				K-1,
				ZAZ,
				m+1,
				IPIV,
				ZAZVT,
				m+1,
				WORK,
				-1);
		WORK = (double *)realloc(WORK, WORK[0]*sizeof(double));
		status = dsysv(
				'L',
				m+1,
				K-1,
				ZAZ,
				m+1,
				IPIV,
				ZAZVT,
				m+1,
				WORK,
				sizeof(WORK)/sizeof(double));
		if (status != 0)
			fprintf(stderr, "Received nonzero status from "
					"dsysv: %i\n", status);
		free(WORK);
		free(IPIV);
	}

	// Return to Row-major order. The matrix ZAZVT contains the solution 
	// after the dposv/dsysv call.
	for (i=0; i<m+1; i++)
		for (j=0; j<K-1; j++)
			ZAZV[i*(K-1)+j] = ZAZVT[j*(m+1)+i];

	// Store the previous V in Vbar, assign the new V
	// (which is stored in ZAZVT) to the model, and give ZAZVT the
	// address of Vbar. This should ensure that we keep
	// re-using assigned memory instead of reallocating at every
	// update.
	/* See this answer: http://stackoverflow.com/q/13246615/
	 * For now we'll just do it by value until the rest is figured out.
	ptr = model->Vbar;
	model->Vbar = model->V;
	model->V = ZAZVT;
	ZAZVT = ptr;
	*/

	for (i=0; i<m+1; i++) {
		for (j=0; j<K-1; j++) {
			matrix_set(
					model->Vbar,
					K-1,
					i,
					j,
					matrix_get(model->V, K-1, i, j));
			matrix_set(
					model->V,
					K-1,
					i,
					j,
					matrix_get(ZAZV, K-1, i, j));
		}
	}
}
