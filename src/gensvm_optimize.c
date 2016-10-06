/**
 * @file gensvm_optimize.c
 * @author Gertjan van den Burg
 * @date August 9, 2013
 * @brief Main functions for training the GenSVM solution.
 *
 * @details
 * Contains update and loss functions used to actually find
 * the optimal V.
 *
 */

#include "gensvm_optimize.h"

/**
 * Maximum number of iterations of the algorithm.
 */
#define MAX_ITER 1000000000

/**
 * Iteration frequency with which to print to stdout
 */
#define PRINT_ITER 100

/**
 * @brief The main training loop for GenSVM
 *
 * @details
 * This function is the main training function. This function
 * handles the optimization of the model with the given model parameters, with
 * the data given. On return the matrix GenModel::V contains the optimal
 * weight matrix.
 *
 * In this function, step doubling is used in the majorization algorithm after
 * a burn-in of 50 iterations.
 *
 * @param[in,out] 	model 	the GenModel to be trained. Contains optimal
 * 				V on exit.
 * @param[in] 		data 	the GenData to train the model with.
 */
void gensvm_optimize(struct GenModel *model, struct GenData *data)
{
	long it = 0;
	double L, Lbar;

	long n = model->n;
	long m = model->m;
	long K = model->K;

	// initialize the workspace
	struct GenWork *work = gensvm_init_work(model);

	// print some info on the dataset and model configuration
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

	// compute necessary simplex vectors
	gensvm_simplex(model);
	gensvm_simplex_diff(model);

	// get initial loss
	L = gensvm_get_loss(model, data, work);
	Lbar = L + 2.0*model->epsilon*L;

	// run main loop
	while ((it < MAX_ITER) && (Lbar - L)/L > model->epsilon)
	{
		// ensures V contains newest V and Vbar contains V from
		// previous
		gensvm_get_update(model, data, work);
		if (it > 50)
			gensvm_step_doubling(model);

		Lbar = L;
		L = gensvm_get_loss(model, data, work);

		if (it%PRINT_ITER == 0)
			note("iter = %li, L = %15.16f, Lbar = %15.16f, "
			     "reldiff = %15.16f\n", it, L, Lbar, (Lbar - L)/L);
		it++;
	}

	// print warnings if necessary
	if (L > Lbar)
		err("[GenSVM Warning]: Negative step occurred in "
				"majorization.\n");
	if (it >= MAX_ITER)
		err("[GenSVM Warning]: maximum number of iterations "
				"reached.\n");

	// print final iteration count and loss
	note("Optimization finished, iter = %li, loss = %15.16f, "
			"rel. diff. = %15.16f\n", it-1, L,
			(Lbar - L)/L);

	// compute and print the number of SVs in the model
	note("Number of support vectors: %li\n", gensvm_num_sv(model));

	// store the training error in the model
	model->training_error = (Lbar - L)/L;

	// free the workspace
	gensvm_free_work(work);
}

/**
 * @brief Calculate the current value of the loss function
 *
 * @details
 * The current loss function value is calculated based on the matrix V in the
 * given model. Note that the matrix ZV is passed explicitly to avoid having
 * to reallocate memory at every step.
 *
 * @param[in] 		model 	GenModel structure which holds the current
 * 				estimate V
 * @param[in]  		data 	GenData structure
 * @param[in] 		work 	allocated workspace with the ZV matrix to use
 * @returns 			the current value of the loss function
 */
double gensvm_get_loss(struct GenModel *model, struct GenData *data, 
		struct GenWork *work)
{
	long i, j;
	long n = model->n;
	long K = model->K;
	long m = model->m;

	double value, rowvalue, loss = 0.0;

	gensvm_calculate_errors(model, data, work->ZV);
	gensvm_calculate_huber(model);

	for (i=0; i<n; i++) {
		rowvalue = 0;
		value = 0;
		for (j=0; j<K; j++) {
			if (j == (data->y[i]-1))
				continue;
			value = matrix_get(model->H, K, i, j);
			value = pow(value, model->p);
			rowvalue += value;
		}
		rowvalue = pow(rowvalue, 1.0/(model->p));
		rowvalue *= model->rho[i];
		loss += rowvalue;
	}
	loss /= ((double) n);

	value = 0;
	for (i=1; i<m+1; i++) {
		for (j=0; j<K-1; j++) {
			value += pow(matrix_get(model->V, K-1, i, j), 2.0);
		}
	}
	loss += model->lambda * value;

	return loss;
}

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
		h = matrix_get(model->H, model->K, i, j);
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
		h = matrix_get(model->H, model->K, i, j);
		value += (h > 0) ? 1 : 0;
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
	double q = matrix_get(model->Q, model->K, i, j);
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
	double q = matrix_get(model->Q, model->K, i, j);

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
	long j,
	     K = model->K;
	double omega, a, b_aq = 0.0,
	       alpha = 0.0;
	double *uu_row;
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
		uu_row = &model->UU[((data->y[i]-1)*K+j)*(K-1)];
		cblas_daxpy(K-1, b_aq, uu_row, 1, beta, 1);

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
	double alpha;

	long n = model->n;
	long m = model->m;
	long K = model->K;

	gensvm_reset_work(work);

	// generate Z'*A*Z and Z'*B by rank 1 operations
	for (i=0; i<n; i++) {
		alpha = gensvm_get_alpha_beta(model, data, i, work->beta);

		// calculate row of matrix LZ, which is a scalar 
		// multiplication of sqrt(alpha_i) and row z_i' of Z
		cblas_daxpy(m+1, sqrt(alpha), &data->Z[i*(m+1)], 1, 
				&work->LZ[i*(m+1)], 1);

		// rank 1 update of matrix Z'*B
		// Note: LDA is the second dimension of ZB because of
		// Row-Major order
		cblas_dger(CblasRowMajor, m+1, K-1, 1, &data->Z[i*(m+1)], 1,
				work->beta, 1, work->ZB, K-1);
	}

	// calculate Z'*A*Z by symmetric multiplication of LZ with itself 
	// (ZAZ = (LZ)' * (LZ)
	cblas_dsyrk(CblasRowMajor, CblasUpper, CblasTrans, m+1, n, 1.0,
			work->LZ, m+1, 0.0, work->ZAZ, m+1);

	// Calculate right-hand side of system we want to solve
	// dsymm performs ZB := 1.0 * (ZAZ) * Vbar + 1.0 * ZB
	// the right-hand side is thus stored in ZB after this call
	// Note: LDB and LDC are second dimensions of the matrices due to
	// Row-Major order
	cblas_dsymm(CblasRowMajor, CblasLeft, CblasUpper, m+1, K-1, 1, 
			work->ZAZ, m+1, model->V, K-1, 1.0, work->ZB, K-1);

	// Calculate left-hand side of system we want to solve
	// Add lambda to all diagonal elements except the first one. Recall
	// that ZAZ is of size m+1 and is symmetric.
	for (i=m+2; i<=m*(m+2); i+=m+2)
		work->ZAZ[i] += model->lambda;

	// Lapack uses column-major order, so we transform the ZB matrix to
	// correspond to this.
	for (i=0; i<m+1; i++)
		for (j=0; j<K-1; j++)
			work->ZBc[j*(m+1)+i] = work->ZB[i*(K-1)+j];

	// Solve the system using dposv. Note that above the upper triangular
	// part has always been used in row-major order for ZAZ. This
	// corresponds to the lower triangular part in column-major order.
	status = dposv('L', m+1, K-1, work->ZAZ, m+1, work->ZBc, m+1);

	// Use dsysv as fallback, for when the ZAZ matrix is not positive
	// semi-definite for some reason (perhaps due to rounding errors).
	// This step shouldn't be necessary but is included for safety.
	if (status != 0) {
		err("[GenSVM Warning]: Received nonzero status from "
				"dposv: %i\n", status);
		int *IPIV = Malloc(int, m+1);
		double *WORK = Malloc(double, 1);
		status = dsysv('L', m+1, K-1, work->ZAZ, m+1, IPIV, work->ZBc, 
				m+1, WORK, -1);

		int LWORK = WORK[0];
		WORK = Realloc(WORK, double, LWORK);
		status = dsysv('L', m+1, K-1, work->ZAZ, m+1, IPIV, work->ZBc, 
				m+1, WORK, LWORK);
		if (status != 0)
			err("[GenSVM Warning]: Received nonzero "
					"status from dsysv: %i\n", status);
		free(WORK);
		free(IPIV);
	}

	// the solution is now stored in ZBc, in column-major order. Here we
	// convert this back to row-major order
	for (i=0; i<m+1; i++)
		for (j=0; j<K-1; j++)
			work->ZB[i*(K-1)+j] = work->ZBc[j*(m+1)+i];

	// copy the old V to Vbar and the new solution to V
	for (i=0; i<m+1; i++) {
		for (j=0; j<K-1; j++) {
			matrix_set(model->Vbar, K-1, i, j,
					matrix_get(model->V, K-1, i, j));
			matrix_set(model->V, K-1, i, j,
					matrix_get(work->ZB, K-1, i, j));
		}
	}
}

/**
 * @brief Use step doubling
 *
 * @details
 * Step doubling can be used to speed up the maorization algorithm. Instead of
 * using the value at the minimimum of the majorization function, the value
 * ``opposite'' the majorization point is used. This can essentially cut the
 * number of iterations necessary to reach the minimum in half.
 *
 * @param[in] 	model	GenModel containing the augmented parameters
 */
void gensvm_step_doubling(struct GenModel *model)
{
	long i, j;
	double value;

	long m = model->m;
	long K = model->K;

	for (i=0; i<m+1; i++) {
		for (j=0; j<K-1; j++) {
			matrix_mul(model->V, K-1, i, j, 2.0);
			value = - matrix_get(model->Vbar, K-1, i, j);
			matrix_add(model->V, K-1, i, j, value);
		}
	}
}

/**
 * @brief Calculate the Huber hinge errors
 *
 * @details
 * For each of the scalar errors in Q the Huber hinge errors are
 * calculated. The Huber hinge is here defined as
 * @f[
 * 	h(q) =
 * 		\begin{dcases}
 * 			1 - q - \frac{\kappa + 1}{2} & \text{if } q \leq
 * 			-\kappa \\
 * 			\frac{1}{2(\kappa + 1)} ( 1 - q)^2 & \text{if } q \in
 * 			(-\kappa, 1] \\
 * 			0 & \text{if } q > 1
 * 		\end{dcases}
 * @f]
 *
 * @param[in,out] model 	the corresponding GenModel
 */
void gensvm_calculate_huber(struct GenModel *model)
{
	long i, j;
	double q, value;

	for (i=0; i<model->n; i++) {
		for (j=0; j<model->K; j++) {
			q = matrix_get(model->Q, model->K, i, j);
			value = 0.0;
			if (q <= -model->kappa) {
				value = 1.0 - q - (model->kappa+1.0)/2.0;
			} else if (q <= 1.0) {
				value = 1.0/(2.0*model->kappa+2.0)*pow(1.0 - q,
					       	2.0);
			}
			matrix_set(model->H, model->K, i, j, value);
		}
	}
}

/**
 * @brief Calculate the scalar errors
 *
 * @details
 * Calculate the scalar errors q based on the current estimate of V, and
 * store these in Q. It is assumed that the memory for Q has already been
 * allocated. In addition, the matrix ZV is calculated here. It is assigned
 * to a pre-allocated block of memory, which is passed to this function.
 *
 * @param[in,out] 	model 	the corresponding GenModel
 * @param[in] 		data 	the corresponding GenData
 * @param[in,out] 	ZV 	a pointer to a memory block for ZV. On exit
 * 				this block is updated with the new ZV matrix
 * 				calculated with GenModel::V
 */
void gensvm_calculate_errors(struct GenModel *model, struct GenData *data,
		double *ZV)
{
	long i, j;
	double q, *uu_row;

	long n = model->n;
	long m = model->m;
	long K = model->K;

	cblas_dgemm(
			CblasRowMajor,
			CblasNoTrans,
			CblasNoTrans,
			n,
			K-1,
			m+1,
			1.0,
			data->Z,
			m+1,
			model->V,
			K-1,
			0,
			ZV,
			K-1);

	for (i=0; i<n; i++) {
		for (j=0; j<K; j++) {
			if (j == (data->y[i]-1))
				continue;
			uu_row = &model->UU[((data->y[i]-1)*K+j)*(K-1)];
			q = cblas_ddot(K-1, &ZV[i*(K-1)], 1, uu_row, 1);
			matrix_set(model->Q, K, i, j, q);
		}
	}
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
