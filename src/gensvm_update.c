/**
 * @file gensvm_update.c
 * @author G.J.J. van den Burg
 * @date 2016-10-14
 * @brief Functions for getting an update of the majorization algorithm
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

#include "gensvm_update.h"

/**
 * Number of rows in a single block for the ZAZ calculation in 
 * gensvm_get_ZAZ_ZB_sparse().
 */
#ifndef GENSVM_BLOCK_SIZE
  #define GENSVM_BLOCK_SIZE 512
#endif

/**
 * @brief Calculate the value of omega for a single instance
 *
 * @details
 * This function calculates the value of the @f$ \omega_i @f$ variable for a
 * single instance, where
 * @f[
 * 	\omega_i = \frac{1}{p} \left( \sum_{j \neq y_i} h^p\left(
 * 	\overline{q}_i^{(y_i j)} \right)  \right)^{1/p-1}
 * @f]
 * Note that his function uses the precalculated values from GenModel::H and
 * GenModel::R to speed up the computation.
 *
 * @param[in] 	model 	GenModel structure with the current model
 * @param[in] 	data 	GenData structure with the data (used for y)
 * @param[in] 	i 	index of the instance for which to calculate omega
 * @returns 		the value of omega for instance i
 *
 */
double gensvm_calculate_omega(struct GenModel *model, struct GenData *data,
		long i)
{
	long j;
	double h, omega = 0.0,
	       p = model->p;

	for (j=0; j<model->K; j++) {
		if (j == (data->y[i]-1))
			continue;
		h = matrix_get(model->H, model->n, model->K, i, j);
		omega += pow(h, p);
	}
	omega = (1.0/p)*pow(omega, 1.0/p - 1.0);

	return omega;
}

/**
 * @brief Check if we can do simple majorization for a given instance
 *
 * @details
 * A simple majorization is possible if at most one of the Huberized hinge
 * errors is nonzero for an instance. This is checked here. For this we
 * compute the product of the Huberized error for all @f$j \neq y_i@f$ and
 * check if strictly less than 2 are nonzero. See also the @ref update_math.
 *
 * @param[in] 	model 	GenModel structure with the current model
 * @param[in] 	data 	GenData structure with the data (used for y)
 * @param[in] 	i 	index of the instance for which to check
 * @returns 		whether or not we can do simple majorization
 *
 */
bool gensvm_majorize_is_simple(struct GenModel *model, struct GenData *data,
		long i)
{
	long j;
	double h, value = 0;
	for (j=0; j<model->K; j++) {
		if (j == (data->y[i]-1))
			continue;
		h = matrix_get(model->H, model->n, model->K, i, j);
		value += h > 0;
		if (value > 1)
			return false;
	}
	return true;
}

/**
 * @brief Compute majorization coefficients for non-simple instance
 *
 * @details
 * In this function we compute the majorization coefficients needed for an
 * instance with a non-simple majorization (@f$\varepsilon_i = 0@f$). In this
 * function, we distinguish a number of cases depending on the value of
 * GenModel::p and the respective value of @f$\overline{q}_i^{(y_ij)}@f$. Note
 * that the linear coefficient is of the form @f$b - a\overline{q}@f$, but
 * often the second term is included in the definition of @f$b@f$, so it can
 * be optimized out. The output argument \p b_aq contains this difference
 * therefore in one go. More details on this function can be found in the @ref
 * update_math. See also gensvm_calculate_ab_simple().
 *
 * @param[in] 	model 	GenModel structure with the current model
 * @param[in] 	i 	index for the instance
 * @param[in] 	j 	index for the class
 * @param[out] 	*a 	output argument for the quadratic coefficient
 * @param[out]  *b_aq 	output argument for the linear coefficient.
 *
 */
void gensvm_calculate_ab_non_simple(struct GenModel *model, long i, long j,
		double *a, double *b_aq)
{
	double q = matrix_get(model->Q, model->n, model->K, i, j);
	double p = model->p;
	double kappa = model->kappa;
	const double a2g2 = 0.25*p*(2.0*p - 1.0)*pow((kappa+1.0)/2.0,p-2.0);

	if (2.0 - model->p < 1e-2) {
		if (q <= - kappa) {
			*b_aq = 0.5 - kappa/2.0 - q;
		} else if ( q <= 1.0) {
			*b_aq = pow(1.0 - q, 3.0)/(2.0*pow(kappa + 1.0, 2.0));
		} else {
			*b_aq = 0;
		}
		*a = 1.5;
	} else {
		if (q <= (p + kappa - 1.0)/(p - 2.0)) {
			*a = 0.25*pow(p, 2.0)*pow(0.5 - kappa/2.0 - q, p - 2.0);
		} else if (q <= 1.0) {
			*a = a2g2;
		} else {
			*a = 0.25*pow(p, 2.0)*pow((p/(p - 2.0))*(0.5 -
						kappa/2.0 - q), p - 2.0);
			*b_aq = (*a)*(2.0*q + kappa - 1.0)/(p - 2.0) +
				0.5*p*pow(p/(p - 2.0)*(0.5 - kappa/2.0 - q),
						p - 1.0);
		}
		if (q <= -kappa) {
			*b_aq = 0.5*p*pow(0.5 - kappa/2.0 - q, p - 1.0);
		} else if ( q <= 1.0) {
			*b_aq = p*pow(1.0 - q, 2.0*p - 1.0)/pow(2*kappa+2.0, p);
		}
	}
}

/**
 * @brief Compute majorization coefficients for simple instances
 *
 * @details
 * In this function we compute the majorization coefficients needed for an
 * instance with a simple majorization. This corresponds to the non-simple
 * majorization for the case where GenModel::p equals 1. Due to this condition
 * the majorization coefficients are quite simple to compute.  Note that the
 * linear coefficient of the majorization is of the form @f$b -
 * a\overline{q}@f$, but often the second term is included in the definition
 * of @f$b@f$, so it can be optimized out. For more details see the @ref
 * update_math, and gensvm_calculate_ab_non_simple().
 *
 * @param[in] 	model 	GenModel structure with the current model
 * @param[in] 	i 	index for the instance
 * @param[in] 	j 	index for the class
 * @param[out] 	*a 	output argument for the quadratic coefficient
 * @param[out] 	*b_aq 	output argument for the linear coefficient
 *
 */
void gensvm_calculate_ab_simple(struct GenModel *model, long i, long j,
		double *a, double *b_aq)
{
	double q = matrix_get(model->Q, model->n, model->K, i, j);

	if (q <= - model->kappa) {
		*a = 0.25/(0.5 - model->kappa/2.0 - q);
		*b_aq = 0.5;
	} else if (q <= 1.0) {
		*a = 1.0/(2.0*model->kappa + 2.0);
		*b_aq = (1.0 - q)*(*a);
	} else {
		*a = -0.25/(0.5 - model->kappa/2.0 - q);
		*b_aq = 0;
	}
}

/**
 * @brief Compute the alpha_i and beta_i for an instance
 *
 * @details
 * This computes the @f$\alpha_i@f$ value for an instance, and simultaneously
 * updating the row of the B matrix corresponding to that
 * instance (the @f$\boldsymbol{\beta}_i'@f$). The advantage of doing this at
 * the same time is that we can compute the a and b values simultaneously in
 * the gensvm_calculate_ab_simple() and gensvm_calculate_ab_non_simple()
 * functions.
 *
 * The computation is done by first checking whether simple majorization is
 * possible for this instance. If so, the @f$\omega_i@f$ value is set to 1.0,
 * otherwise this value is computed. If simple majorization is possible, the
 * coefficients a and b_aq are computed by gensvm_calculate_ab_simple(),
 * otherwise they're computed by gensvm_calculate_ab_non_simple(). Next, the
 * beta_i updated through the efficient BLAS daxpy function, and part of the
 * value of @f$\alpha_i@f$ is computed. The final value of @f$\alpha_i@f$ is
 * returned.
 *
 * @param[in] 		model 	GenModel structure with the current model
 * @param[in] 		data 	GenData structure with the data
 * @param[in] 		i 	index of the instance to update
 * @param[out] 		beta	beta vector of linear coefficients (assumed to
 * 				be allocated elsewhere, initialized here)
 * @returns 			the @f$\alpha_i@f$ value of this instance
 *
 */
double gensvm_get_alpha_beta(struct GenModel *model, struct GenData *data,
		long i, double *beta)
{
	bool simple;
	long j, K = model->K;
	double omega, a, b_aq = 0.0,
	       alpha = 0.0;
	double *uu_row = NULL;
	const double in = 1.0/((double) model->n);

	simple = gensvm_majorize_is_simple(model, data, i);
	omega = simple ? 1.0 : gensvm_calculate_omega(model, data, i);

	Memset(beta, double, K-1);
	for (j=0; j<K; j++) {
		// skip the class y_i = k
		if (j == (data->y[i]-1))
			continue;

		// calculate the a_ijk and (b_ijk - a_ijk q_i^(kj)) values
		if (simple) {
			gensvm_calculate_ab_simple(model, i, j, &a, &b_aq);
		} else {
			gensvm_calculate_ab_non_simple(model, i, j, &a, &b_aq);
		}

		// daxpy on beta and UU
		// daxpy does: y = a*x + y
		// so y = beta, UU_row = x, a = factor
		b_aq *= model->rho[i] * omega * in;
		uu_row = &matrix_get(model->UU, K*K, K-1, 
				(data->y[i]-1)*K+j, 0);
		#if MAJOR_ORDER == 'r'
		cblas_daxpy(K-1, b_aq, uu_row, 1, beta, 1);
		#else
		cblas_daxpy(K-1, b_aq, uu_row, K*K, beta, 1);
		#endif

		// increment Avalue
		alpha += a;
	}
	alpha *= omega * model->rho[i] * in;
	return alpha;
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
 * Because the function gensvm_get_update() is always called after a call to
 * gensvm_get_loss() with the same GenModel::V, it is unnecessary to calculate
 * the updated errors GenModel::Q and GenModel::H here too. This saves on
 * computation time.
 *
 * In calculating the majorization coefficients we calculate the elements of a
 * diagonal matrix A with elements
 * @f[
 * 	A_{i, i} = \frac{1}{n} \rho_i \sum_{j \neq k} \left[
 * 		\varepsilon_i a_{ijk}^{(1)} + (1 - \varepsilon_i) \omega_i
 * 		a_{ijk}^{(p)} \right],
 * @f]
 * where @f$ k = y_i @f$.
 * Since this matrix is only used to calculate the matrix @f$ Z' A Z @f$, it
 * is efficient to update a matrix ZAZ through consecutive rank 1 updates with
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
 * Consider using CblasColMajor everywhere
 *
 * @param[in,out] 	model 	model to be updated
 * @param[in] 		data 	data used in model
 * @param[in] 		work 	allocated workspace to use
 */
void gensvm_get_update(struct GenModel *model, struct GenData *data,
		struct GenWork *work)
{
	int status;
	long i, j;

	long m = model->m;
	long K = model->K;

	// compute the ZAZ and ZB matrices
	gensvm_get_ZAZ_ZB(model, data, work);

	// Calculate right-hand side of system we want to solve
	// dsymm performs ZB := 1.0 * (ZAZ) * Vbar + 1.0 * ZB
	// the right-hand side is thus stored in ZB after this call

	#if MAJOR_ORDER == 'r'
	// Note: LDB and LDC are second dimensions of the matrices due to
	// Row-Major order
	cblas_dsymm(CblasRowMajor, CblasLeft, CblasUpper, m+1, K-1, 1,
			work->ZAZ, m+1, model->V, K-1, 1.0, work->ZB, K-1);
	#else
	cblas_dsymm(CblasColMajor, CblasLeft, CblasUpper, m+1, K-1, 1.0,
			work->ZAZ, m+1, model->V, m+1, 1.0, work->ZB, m+1);
	#endif

	// Calculate left-hand side of system we want to solve
	// Add lambda to all diagonal elements except the first one. Recall
	// that ZAZ is of size m+1 and is symmetric.
	for (i=m+2; i<=m*(m+2); i+=m+2)
		work->ZAZ[i] += model->lambda;

	#if MAJOR_ORDER == 'r'
	// Lapack uses column-major order, so we transform the ZB matrix to
	// correspond to this.
	for (i=0; i<m+1; i++)
		for (j=0; j<K-1; j++)
			work->ZBc[j*(m+1)+i] = work->ZB[i*(K-1)+j];
	#endif

	// Solve the system using dposv.
	#if MAJOR_ORDER == 'r'
	// Note that above the upper triangular part has always been used in 
	// row-major order for ZAZ. This corresponds to the lower triangular 
	// part in column-major order.
	status = dposv('L', m+1, K-1, work->ZAZ, m+1, work->ZBc, m+1);
	#else
	status = dposv('U', m+1, K-1, work->ZAZ, m+1, work->ZB, m+1);
	#endif

	// Use dsysv as fallback, for when the ZAZ matrix is not positive
	// semi-definite for some reason (perhaps due to rounding errors).
	// This step shouldn't be necessary but is included for safety.
	if (status != 0) {
		gensvm_error("[GenSVM Warning]: Received nonzero status from "
				"dposv: %i\n", status);
		int *IPIV = Malloc(int, m+1);
		double *WORK = Malloc(double, 1);
		#if MAJOR_ORDER == 'r'
		status = dsysv('L', m+1, K-1, work->ZAZ, m+1, IPIV, work->ZBc,
				m+1, WORK, -1);
		#else
		status = dsysv('U', m+1, K-1, work->ZAZ, m+1, IPIV, work->ZB,
				m+1, WORK, -1);
		#endif

		int LWORK = WORK[0];
		WORK = Realloc(WORK, double, LWORK);

		#if MAJOR_ORDER == 'r'
		status = dsysv('L', m+1, K-1, work->ZAZ, m+1, IPIV, work->ZBc,
				m+1, WORK, LWORK);
		#else
		status = dsysv('U', m+1, K-1, work->ZAZ, m+1, IPIV, work->ZB,
				m+1, WORK, LWORK);
		#endif
		if (status != 0)
			gensvm_error("[GenSVM Warning]: Received nonzero "
					"status from dsysv: %i\n", status);

		free(WORK);
		WORK = NULL;
		free(IPIV);
		IPIV = NULL;
	}

	#if MAJOR_ORDER == 'r'
	// the solution is now stored in ZBc, in column-major order. Here we
	// convert this back to row-major order
	for (i=0; i<m+1; i++)
		for (j=0; j<K-1; j++)
			work->ZB[i*(K-1)+j] = work->ZBc[j*(m+1)+i];
	#endif

	// copy the old V to Vbar and the new solution to V
	for (i=0; i<m+1; i++) {
		for (j=0; j<K-1; j++) {
			matrix_set(model->Vbar, m+1, K-1, i, j,
					matrix_get(model->V, m+1, K-1, i, j));
			matrix_set(model->V, m+1, K-1, i, j,
					matrix_get(work->ZB, m+1, K-1, i, j));
		}
	}
}

/**
 * @brief Calculate Z'*A*Z and Z'*B for dense matrices
 *
 * @details
 * This function calculates the matrices Z'*A*Z and Z'*B for the case where Z
 * is stored as a dense matrix. It calculates the Z'*A*Z product by
 * constructing a matrix LZ = (A^(1/2) * Z), and calculating (LZ)'*(LZ) with
 * the BLAS dsyrk function. The matrix Z'*B is calculated with successive
 * rank-1 updates using the BLAS dger function. These functions came out as
 * the most efficient way to do these computations in several simulation
 * studies.
 *
 * @param[in] 		model 	a GenModel holding the current model
 * @param[in] 		data 	a GenData with the data
 * @param[in,out] 	work 	an allocated GenWork structure, contains
 * 				updated ZAZ and ZB matrices on exit.
 */
void gensvm_get_ZAZ_ZB_dense(struct GenModel *model, struct GenData *data,
		struct GenWork *work)
{
	long i;
	double alpha, sqalpha;

	long n = model->n;
	long m = model->m;
	long K = model->K;

	// generate Z'*A*Z and Z'*B by rank 1 operations
	for (i=0; i<n; i++) {
		alpha = gensvm_get_alpha_beta(model, data, i, work->beta);

		// calculate row of matrix LZ, which is a scalar
		// multiplication of sqrt(alpha_i) and row z_i' of Z
		// Note that we use the fact that the first column of Z is
		// always 1, by only computing the product for m values and
		// copying the first element over.
		sqalpha = sqrt(alpha);

		#if MAJOR_ORDER == 'r'
		work->LZ[i*(m+1)] = sqalpha;
		cblas_daxpy(m, sqalpha, &data->Z[i*(m+1)+1], 1,
				&work->LZ[i*(m+1)+1], 1);
		#else
		work->LZ[i] = sqalpha;
		cblas_daxpy(m, sqalpha, &data->Z[i+n], n, &work->LZ[i+n], n);
		#endif

		// rank 1 update of matrix Z'*B
		#if MAJOR_ORDER == 'r'
		// Note: LDA is the second dimension of ZB because of
		// Row-Major order
		cblas_dger(CblasRowMajor, m+1, K-1, 1, &data->Z[i*(m+1)], 1,
				work->beta, 1, work->ZB, K-1);
		#else
		cblas_dger(CblasColMajor, m+1, K-1, 1, &data->Z[i], n,
				work->beta, 1, work->ZB, m+1);
		#endif

	}

	// calculate Z'*A*Z by symmetric multiplication of LZ with itself
	// (ZAZ = (LZ)' * (LZ)
	#if MAJOR_ORDER == 'r'
	cblas_dsyrk(CblasRowMajor, CblasUpper, CblasTrans, m+1, n, 1.0,
			work->LZ, m+1, 0.0, work->ZAZ, m+1);
	#else
	cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, m+1, n, 1.0,
			work->LZ, n, 0.0, work->ZAZ, m+1);
	#endif
}

void gensvm_get_ZAZ_ZB_sparse(struct GenModel *model, struct GenData *data,
		struct GenWork *work)
{
	if (data->spZ->type == CSR)
		gensvm_get_ZAZ_ZB_sparse_csr(model, data, work);
	else
		gensvm_get_ZAZ_ZB_sparse_csc(model, data, work);
}


/**
 * @brief Calculate Z'*A*Z and Z'*B for sparse matrices
 *
 * @details
 * This function calculates the matrices Z'*A*Z and Z'*B for the case where Z 
 * is stored as a CSR sparse matrix (GenSparse structure). It computes only 
 * the products of the Z'*A*Z matrix that need to be computed, and updates the 
 * Z'*B matrix row-wise for each non-zero element of a row of Z, using a BLAS 
 * daxpy call.
 *
 * This function calculates the matrix product Z'*A*Z in separate blocks, 
 * based on the number of rows defined in the GENSVM_BLOCK_SIZE variable. This 
 * is done to improve numerical precision for very large datasets. Due to 
 * rounding errors, precision can become an issue for these large datasets, 
 * when separate blocks are used and added to the result separately, this can 
 * be alleviated a little bit. See also: http://stackoverflow.com/q/40286989
 *
 * @sa
 * gensvm_get_ZAZ_ZB()
 * gensvm_get_ZAZ_ZB_dense()
 *
 * @param[in] 		model 	a GenModel holding the current model
 * @param[in] 		data 	a GenData with the data
 * @param[in,out] 	work 	an allocated GenWork structure, contains
 * 				updated ZAZ and ZB matrices on exit.
 */
void gensvm_get_ZAZ_ZB_sparse_csr(struct GenModel *model, struct GenData *data,
		struct GenWork *work)
{
	long *Zi = NULL,
	     *Zj = NULL;
	long b, i, j, k, K, kk, b_start, b_end, blk, blk_start, blk_end,
	     rem_size, n_blocks, n_row = data->spZ->n_row,
	     n_col = data->spZ->n_col;
	double temp, alpha, z_ij, *vals = NULL;

	K = model->K;
	Zi = data->spZ->ix;
	Zj = data->spZ->jx;
	vals = data->spZ->values;

	// calculate ZAZ using blocks of rows of Z. This helps avoiding 
	// rounding errors, which increases precision, and in turn helps 
	// convergence of the IM algorithm.
	// see also: http://stackoverflow.com/q/40286989/
	n_blocks = floor(n_row / GENSVM_BLOCK_SIZE);
	rem_size = n_row % GENSVM_BLOCK_SIZE;

	for (blk=0; blk<=n_blocks; blk++) {
		blk_start = blk * GENSVM_BLOCK_SIZE;
		blk_end = blk_start;
		blk_end += (blk == n_blocks) ? rem_size : GENSVM_BLOCK_SIZE;

		Memset(work->tmpZAZ, double, n_col*n_col);
		for (i=blk_start; i<blk_end; i++) {
			alpha = gensvm_get_alpha_beta(model, data, i, 
					work->beta);
			b_start = Zi[i];
			b_end = Zi[i+1];

			for (b=b_start; b<b_end; b++) {
				j = Zj[b];
				z_ij = vals[b];
				cblas_daxpy(K-1, z_ij, work->beta, 1,
						&work->ZB[j*(K-1)], 1);

				z_ij *= alpha;
				for (kk=b; kk<b_end; kk++) {
					matrix_add(work->tmpZAZ, n_row, 
							n_col, j, Zj[kk], 
							z_ij*vals[kk]);
				}
			}
		}

		// copy the intermediate results over to the actual ZAZ matrix
		for (j=0; j<n_col; j++) {
			for (k=j; k<n_col; k++) {
				temp = matrix_get(work->tmpZAZ, n_col, n_col, j, k);
				matrix_add(work->ZAZ, n_col, n_col, j, k, temp);
			}
		}
	}
}


void gensvm_get_ZAZ_ZB_sparse_csc(struct GenModel *model, struct GenData *data,
		struct GenWork *work)
{
	long i, j, k, l, aa, aa_start, aa_end, bb, bb_start, bb_end;
	long *Zi = NULL,
	     *Zj = NULL;
	double z_ik, z_jl, alpha, b_il;
	double *vals = NULL;
	double *beta = NULL;

	long n = model->n;
	long m = model->m;
	long K = model->K;

	Zi = data->spZ->ix;
	Zj = data->spZ->jx;

	beta = Malloc(double, K-1);
	vals = data->spZ->values;

	// compute A and B
	for (i=0; i<n; i++) {
		work->A[i] = gensvm_get_alpha_beta(model, data, i, beta);
		for (j=0; j<K-1; j++)
			matrix_set(work->B, n, K-1, i, j, beta[j]);
	}

	// Compute Z'*A*Z and Z'*B
	for (k=0; k<data->spZ->n_col; k++) {
		// note that answer is symmetric
		aa_start = Zi[k];
		aa_end = Zi[k+1];

		for (aa=aa_start; aa<aa_end; aa++) {
			i = Zj[aa];
			z_ik = vals[aa];
			alpha = work->A[i];

			for (l=k; l<data->spZ->n_col; l++) {
				bb_start = Zi[l];
				bb_end = Zi[l+1];

				for (bb=bb_start; bb<bb_end; bb++) {
					j = Zj[bb];
					if (i != j)
						continue;

					z_jl = vals[bb];
					matrix_add(
							work->ZAZ,
							m+1,
							m+1,
							k,
							l,
							alpha * z_ik * z_jl
						    );
				}
			}

			// this part is for Z'*B
			for (l=0; l<K-1; l++) {
				b_il = matrix_get(work->B, n, K-1, i, l);
				matrix_add(work->ZB, m+1, K-1, k, l, 
						z_ik * b_il);
			}
		}
	}

	free(beta);
}


/**
 * @brief Wrapper around calculation of Z'*A*Z and Z'*B for sparse and dense
 *
 * @details
 * This is a wrapper around gensvm_get_ZAZ_ZB_dense() and 
 * gensvm_get_ZAZ_ZB_sparse(). See the documentation of those functions for 
 * more info.
 *
 * @param[in]	 model 	a GenModel struct
 * @param[in]	 data 	a GenData struct
 * @param[in]	 work 	a GenWork struct
 *
 */
void gensvm_get_ZAZ_ZB(struct GenModel *model, struct GenData *data,
		struct GenWork *work)
{
	gensvm_reset_work(work);

	if (data->Z == NULL)
		gensvm_get_ZAZ_ZB_sparse(model, data, work);
	else
		gensvm_get_ZAZ_ZB_dense(model, data, work);
}

/**
 * @brief Solve AX = B where A is symmetric positive definite.
 *
 * @details
 * Solve a linear system of equations AX = B where A is symmetric positive
 * definite. This function is a wrapper for the external  LAPACK routine
 * dposv.
 *
 * @param[in] 		UPLO 	which triangle of A is stored
 * @param[in] 		N 	order of A
 * @param[in] 		NRHS 	number of columns of B
 * @param[in,out] 	A 	double precision array of size (LDA, N). On
 * 				exit contains the upper or lower factor of the
 * 				Cholesky factorization of A.
 * @param[in] 		LDA 	leading dimension of A
 * @param[in,out] 	B 	double precision array of size (LDB, NRHS). On
 * 				exit contains the N-by-NRHS solution matrix X.
 * @param[in] 		LDB 	the leading dimension of B
 * @returns 			info parameter which contains the status of the
 * 				computation:
 * 					- =0: 	success
 * 					- <0: 	if -i, the i-th argument had
 * 						an illegal value
 * 					- >0: 	if i, the leading minor of A
 * 						was not positive definite
 *
 * See the LAPACK documentation at:
 * http://www.netlib.org/lapack/explore-html/dc/de9/group__double_p_osolve.html
 */
int dposv(char UPLO, int N, int NRHS, double *A, int LDA, double *B,
		int LDB)
{
	extern void dposv_(char *UPLO, int *Np, int *NRHSp, double *A,
			int *LDAp, double *B, int *LDBp, int *INFOp);
	int INFO;
	dposv_(&UPLO, &N, &NRHS, A, &LDA, B, &LDB, &INFO);
	return INFO;
}

/**
 * @brief Solve a system of equations AX = B where A is symmetric.
 *
 * @details
 * Solve a linear system of equations AX = B where A is symmetric. This
 * function is a wrapper for the external LAPACK routine dsysv.
 *
 * @param[in] 		UPLO 	which triangle of A is stored
 * @param[in] 		N 	order of A
 * @param[in] 		NRHS 	number of columns of B
 * @param[in,out] 	A 	double precision array of size (LDA, N). On
 * 				exit contains the block diagonal matrix D and
 * 				the multipliers used to obtain the factor U or
 * 				L from the factorization A = U*D*U**T or
 * 				A = L*D*L**T.
 * @param[in] 		LDA 	leading dimension of A
 * @param[in] 		IPIV 	integer array containing the details of D
 * @param[in,out] 	B 	double precision array of size (LDB, NRHS). On
 * 				exit contains the N-by-NRHS matrix X
 * @param[in] 		LDB 	leading dimension of B
 * @param[out] 		WORK 	double precision array of size max(1,LWORK). On
 * 				exit, WORK(1) contains the optimal LWORK
 * @param[in] 		LWORK 	the length of WORK, can be used for determining
 * 				the optimal blocksize for dsystrf.
 * @returns 			info parameter which contains the status of the
 * 				computation:
 * 					- =0: 	success
 * 					- <0: 	if -i, the i-th argument had an
 * 						illegal value
 * 					- >0: 	if i, D(i, i) is exactly zero,
 * 						no solution can be computed.
 *
 * See the LAPACK documentation at:
 * http://www.netlib.org/lapack/explore-html/d6/d0e/group__double_s_ysolve.html
 */
int dsysv(char UPLO, int N, int NRHS, double *A, int LDA, int *IPIV,
		double *B, int LDB, double *WORK, int LWORK)
{
	extern void dsysv_(char *UPLO, int *Np, int *NRHSp, double *A,
			int *LDAp, int *IPIV, double *B, int *LDBp,
			double *WORK, int *LWORK, int *INFOp);
	int INFO;
	dsysv_(&UPLO, &N, &NRHS, A, &LDA, IPIV, B, &LDB, WORK, &LWORK, &INFO);
	return INFO;
}
