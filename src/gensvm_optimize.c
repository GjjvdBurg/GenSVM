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
 * @brief The main training loop for GenSVM
 *
 * @details
 * This function is the main training function. This function
 * handles the optimization of the model with the given model parameters, with
 * the data given. On return the matrix GenModel::V contains the optimal
 * weight matrix.
 *
 * In this function, step doubling is used in the majorization algorithm after
 * a burn-in of 50 iterations. If the training is finished, GenModel::t and
 * GenModel::W are extracted from GenModel::V.
 *
 * @param[in,out] 	model 	the GenModel to be trained. Contains optimal
 * 				V on exit.
 * @param[in] 		data 	the GenData to train the model with.
 */
void gensvm_optimize(struct GenModel *model, struct GenData *data)
{
	long i, j, it = 0;
	double L, Lbar, value;

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

	gensvm_simplex(model->K, model->U);
	gensvm_simplex_diff(model, data);
	gensvm_category_matrix(model, data);

	L = gensvm_get_loss(model, data, ZV);
	Lbar = L + 2.0*model->epsilon*L;

	while ((it < MAX_ITER) && (Lbar - L)/L > model->epsilon)
	{
		// ensure V contains newest V and Vbar contains V from
		// previous
		gensvm_get_update(model, data, B, ZAZ, ZAZV, ZAZVT);
		if (it > 50)
			gensvm_step_doubling(model);

		Lbar = L;
		L = gensvm_get_loss(model, data, ZV);

		if (it%100 == 0)
			note("iter = %li, L = %15.16f, Lbar = %15.16f, "
			     "reldiff = %15.16f\n", it, L, Lbar, (Lbar - L)/L);
		it++;
	}
	if (L > Lbar)
		err("[GenSVM Warning]: Negative step occurred in "
				"majorization.\n");
	if (it >= MAX_ITER)
		err("[GenSVM Warning]: maximum number of iterations "
				"reached.\n");

	note("Optimization finished, iter = %li, loss = %15.16f, "
			"rel. diff. = %15.16f\n", it-1, L,
			(Lbar - L)/L);
	note("Number of support vectors: %li\n", gensvm_num_sv(model));

	model->training_error = (Lbar - L)/L;

	for (i=0; i<K-1; i++)
		model->t[i] = matrix_get(model->V, K-1, 0, i);
	for (i=1; i<m+1; i++) {
		for (j=0; j<K-1; j++) {
			value = matrix_get(model->V, K-1, i, j);
			matrix_set(model->W, K-1, i-1, j, value);
		}
	}
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
 * @param[in] 		model 	GenModel structure which holds the current
 * 				estimate V
 * @param[in]  		data 	GenData structure
 * @param[in,out] 	ZV 	pre-allocated matrix ZV which is updated on
 * 				output
 * @returns 			the current value of the loss function
 */
double gensvm_get_loss(struct GenModel *model, struct GenData *data,
		double *ZV)
{
	long i, j;
	long n = model->n;
	long K = model->K;
	long m = model->m;

	double value, rowvalue, loss = 0.0;

	gensvm_calculate_errors(model, data, ZV);
	gensvm_calculate_huber(model);

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
double gensvm_calculate_omega(struct GenModel *model, long i)
{
	long j;
	double h, r, omega = 0.0,
	       p = model->p;

	for (j=0; j<model->K; j++) {
		h = matrix_get(model->H, model->K, i, j);
		r = matrix_get(model->R, model->K, i, j);
		omega += pow(h, p)*r;
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
bool gensvm_majorize_is_simple(struct GenModel *model, long i)
{
	long j;
	double h, r, value = 0;
	for (j=0; j<model->K; j++) {
		h = matrix_get(model->H, model->K, i, j);
		r = matrix_get(model->R, model->K, i, j);
		value += (h*r > 0) ? 1 : 0;
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
 * @brief Update the relevant elements of the B matrix with given coefficients
 *
 * @details
 * This function updates the row of the B matrix for the given instance \p i,
 * with the provided linear coefficient \p b_aq. Each row of the B matrix is a
 * sum of \p b_aq values for class \p j, multiplied with their respective
 * simplex-difference vectors @f$\boldsymbol{\delta_{y_ij}}'@f$. This function
 * adds to the row of B for a single class \p j, it should therefore be
 * repeatedly called for each class to get the correct outcome. Note that the
 * sum should be done over all @f$j \neq y_i@f$, but since the simplex
 * difference vector @f$\boldsymbol{\delta_{jj}}'@f$ is zero, we do not need
 * an extra restriction for this. For more details see @ref update_math.
 *
 * @param[in] 		model 	GenModel structure with the current model
 * @param[in] 		i 	index of the instance
 * @param[in] 		j 	index for the class
 * @param[in] 		b_aq 	linear coefficient update for this i, j
 * @param[in] 		omega 	value of omega to multiply with
 * @param[in,out] 	B 	matrix of linear coefficients
 *
 */
void gensvm_update_B(struct GenModel *model, long i, long j, double b_aq,
		double omega, double *B)
{
	long k;
	double Bvalue;
	const double factor = (model->rho[i]*omega*b_aq)/((double) model->n);

	for (k=0; k<model->K-1; k++) {
		Bvalue = factor * matrix3_get(model->UU, model->K-1, model->K,
				i, k, j);
		matrix_add(B, model->K-1, i, k, Bvalue);
	}
}

/**
 * @brief Update the B matrix for an instance and compute the alpha_i
 *
 * @details
 * This computes the @f$\alpha_i@f$ value for an instance, while
 * simultaneously updating the row of the B matrix corresponding to that
 * instance. The advantage of doing this at the same time is that we can
 * compute the a and b values simultaneously in the
 * gensvm_calculate_ab_simple() and gensvm_calculate_ab_non_simple()
 * functions.
 *
 * The computation is done by first checking whether simple majorization is
 * possible for this instance. If so, the @f$\omega_i@f$ value is set to 1.0,
 * otherwise this value is computed. If simple majorization is possible, the
 * coefficients a and b_aq are computed by gensvm_calculate_ab_simple(),
 * otherwise they're computed by gensvm_calculate_ab_non_simple(). Next, the
 * relevant row of the B matrix is updated, and part of the value of
 * @f$\alpha_i@f$ is computed. The final value of @f$\alpha_i@f$ is returned.
 *
 * @param[in] 		model 	GenModel structure with the current model
 * @param[in] 		i 	index of the instance to update
 * @param[in,out] 	B 	B matrix of linear coefficients
 * @returns 			the @f$\alpha_i@f$ value of this instance
 *
 */
double gensvm_get_Avalue_update_B(struct GenModel *model, long i, double *B)
{
	bool simple;
	long j;
	double omega, a, b_aq, Avalue = 0.0;
	const double in = 1.0/((double) model->n);

	simple = gensvm_majorize_is_simple(model, i);
	omega = simple ? 1.0 : gensvm_calculate_omega(model, i);

	for (j=0; j<model->K; j++) {
		if (simple) {
			gensvm_calculate_ab_simple(model, i, j, &a, &b_aq);
		} else {
			gensvm_calculate_ab_non_simple(model, i, j, &a, &b_aq);
		}
		gensvm_update_B(model, i, j, b_aq, omega, B);
		Avalue += a*matrix_get(model->R, model->K, i, j);
	}
	Avalue *= omega * in * model->rho[i];
	return Avalue;
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
 * Consider allocating IPIV and WORK at a higher level, they probably don't
 * change much during the iterations.
 *
 * @param [in,out] 	model 	model to be updated
 * @param [in] 		data 	data used in model
 * @param [in] 		B 	pre-allocated matrix used for linear
 * 				coefficients
 * @param [in] 		ZAZ 	pre-allocated matrix used in system
 * @param [in] 		ZAZV 	pre-allocated matrix used in system solving
 * @param [in] 		ZAZVT 	pre-allocated matrix used in system solving
 */
void gensvm_get_update(struct GenModel *model, struct GenData *data,
		double *B, double *ZAZ, double *ZAZV, double *ZAZVT)
{
	int status;
	long i, j;
	double Avalue;
	double value;

	long n = model->n;
	long m = model->m;
	long K = model->K;

	// clear matrices
	Memset(B, double, n*(K-1));
	Memset(ZAZ, double, (m+1)*(m+1));

	for (i=0; i<n; i++) {
		Avalue = gensvm_get_Avalue_update_B(model, i, B);

		// Now we calculate the matrix ZAZ. Since this is
		// guaranteed to be symmetric, we only calculate the
		// upper part of the matrix, and then copy this over
		// to the lower part after all calculations are done.
		// Note that the use of dsyr is faster than dspr, even
		// though dspr uses less memory.
		cblas_dsyr(CblasRowMajor, CblasUpper, m+1, Avalue,
				&data->Z[i*(m+1)], 1, ZAZ, m+1);
	}

	// Calculate the right hand side of the system we
	// want to solve.
	cblas_dsymm(CblasRowMajor, CblasLeft, CblasUpper, m+1, K-1, 1.0, ZAZ,
			m+1, model->V, K-1, 0.0, ZAZV, K-1);

	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m+1, K-1, n, 1.0,
			data->Z, m+1, B, K-1, 1.0, ZAZV, K-1);

	// Add lambda to all diagonal elements except the first one. Recall
	// that ZAZ is of size m+1 and is symmetric.
	for (i=m+2; i<=m*(m+2); i+=m+2)
		ZAZ[i] += model->lambda;

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
	status = dposv('L', m+1, K-1, ZAZ, m+1, ZAZVT, m+1);

	if (status != 0) {
		// This step should not be necessary, as the matrix
		// ZAZ is positive semi-definite by definition. It
		// is included for safety.
		err("[GenSVM Warning]: Received nonzero status from "
				"dposv: %i\n", status);
		int *IPIV = Malloc(int, m+1);
		double *WORK = Malloc(double, 1);
		status = dsysv('L', m+1, K-1, ZAZ, m+1, IPIV, ZAZVT, m+1, WORK,
				-1);

		int LWORK = WORK[0];
		WORK = Realloc(WORK, double, LWORK);
		status = dsysv('L', m+1, K-1, ZAZ, m+1, IPIV, ZAZVT, m+1, WORK,
				LWORK);
		if (status != 0)
			err("[GenSVM Warning]: Received nonzero "
					"status from dsysv: %i\n", status);
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
			value = matrix_get(model->V, K-1, i, j);
			matrix_set(model->Vbar, K-1, i, j, value);
			value = matrix_get(ZAZV, K-1, i, j);
			matrix_set(model->V, K-1, i, j, value);
		}
	}
}

/**
 * @brief Generate the category matrix
 *
 * @details
 * Generate the category matrix R. The category matrix has 1's everywhere
 * except at the column corresponding to the label of instance i, there the
 * element is 0.
 *
 * @param[in,out] 	model 		corresponding GenModel
 * @param[in] 		dataset 	corresponding GenData
 *
 */
void gensvm_category_matrix(struct GenModel *model, struct GenData *data)
{
	long i, j;
	long n = model->n;
	long K = model->K;

	for (i=0; i<n; i++) {
		for (j=0; j<K; j++) {
			if (data->y[i] != j+1)
				matrix_set(model->R, K, i, j, 1.0);
			else
				matrix_set(model->R, K, i, j, 0.0);
		}
	}
}

/**
 * @brief Generate the simplex difference matrix
 *
 * @details
 * The simplex difference matrix is a 3D matrix which is constructed
 * as follows. For each instance i, the difference vectors between the row of
 * the simplex matrix corresponding to the class label of instance i and the
 * other rows of the simplex matrix are calculated. These difference vectors
 * are stored in a matrix, which is one horizontal slice of the 3D matrix.
 *
 * We use the indices i, j, k for the three dimensions n, K-1, K of UU. Then
 * the i,j,k -th element of UU is equal to U(y[i]-1, j) - U(k, j).
 *
 * @param[in,out] 	model 	the corresponding GenModel
 * @param[in] 		data 	the corresponding GenData
 *
 */
void gensvm_simplex_diff(struct GenModel *model, struct GenData *data)
{
	long i, j, k;
	double value;

	long n = model->n;
	long K = model->K;

	for (i=0; i<n; i++) {
		for (j=0; j<K-1; j++) {
			for (k=0; k<K; k++) {
				value = matrix_get(model->U, K-1,
						data->y[i]-1, j);
				value -= matrix_get(model->U, K-1, k, j);
				matrix3_set(model->UU, K-1, K, i, j, k, value);
			}
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
 *
 */
void gensvm_calculate_errors(struct GenModel *model, struct GenData *data,
		double *ZV)
{
	long i, j, k;
	double zv, value;

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

	Memset(model->Q, double, n*K);
	for (i=0; i<n; i++) {
		for (j=0; j<K-1; j++) {
			zv = matrix_get(ZV, K-1, i, j);
			for (k=0; k<K; k++) {
				value = zv * matrix3_get(model->UU, K-1, K, i,
						j, k);
				matrix_add(model->Q, K, i, k, value);
			}
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
