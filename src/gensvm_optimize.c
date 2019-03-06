/**
 * @file gensvm_optimize.c
 * @author G.J.J. van den Burg
 * @date 2013-08-09
 * @brief Main functions for training the GenSVM solution.
 *
 * @details
 * Contains update and loss functions used to actually find
 * the optimal V.
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

#include "gensvm_optimize.h"

/**
 * Iteration frequency with which to print to stdout
 */
#ifndef GENSVM_PRINT_ITER
  #define GENSVM_PRINT_ITER 100
#endif

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
	GenTime t_start, t_stop, t_ipt_start, t_ipt_stop;

	gensvm_py_reset_interrupt_hdl();

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

	Timer(t_start);
	Timer(t_ipt_start);

	// get initial loss
	L = gensvm_get_loss(model, data, work);
	Lbar = L + 2.0*model->epsilon*L;

	// run main loop
	while ((it < model->max_iter) && (Lbar - L)/L > model->epsilon)
	{
		// ensures V contains newest V and Vbar contains V from
		// previous
		gensvm_get_update(model, data, work);
		if (it > 50)
			gensvm_step_doubling(model);

		Lbar = L;
		L = gensvm_get_loss(model, data, work);

		if (it % GENSVM_PRINT_ITER == 0)
			note("iter = %li, L = %15.16f, Lbar = %15.16f, "
			     "reldiff = %15.16f\n", it, L, Lbar, (Lbar - L)/L);
		it++;

		Timer(t_ipt_stop);
		if (gensvm_elapsed_time(&t_ipt_start, &t_ipt_stop) > 2) {
			if (gensvm_py_pending_interrupt()) {
				gensvm_error("[GenSVM Warning]: Received "
						"user interrupt. Stopping.\n");
				break;
			}
			Timer(t_ipt_start);
		}
	}

	Timer(t_stop);

	// status == 0 means training was successful
	model->status = 0;

	// print warnings if necessary
	if (L > Lbar) {
		gensvm_error("[GenSVM Warning]: Negative step occurred in "
				"majorization.\n");
		model->status = 1;
	}

	if (it >= model->max_iter) {
		gensvm_error("[GenSVM Warning]: maximum number of iterations "
				"reached.\n");
		model->status = 2;
	}

	if (gensvm_py_pending_interrupt())
		model->status = 3;

	// print final iteration count and loss
	note("\nOptimization finished, iter = %li, loss = %15.16f, "
			"reldiff = %15.16f\n", it, L,
			(Lbar - L)/L);

	// compute and print the number of SVs in the model
	note("Number of support vectors: %li\n", gensvm_num_sv(model));

	// store the training error in the model
	model->training_error = (Lbar - L)/L;

	// store the iteration count in the model
	model->elapsed_iter = it;

	// store the elapsed time in the model
	model->elapsed_time = gensvm_elapsed_time(&t_start, &t_stop);
	note("Training time: %f\n", model->elapsed_time);

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
			value = matrix_get(model->H, n, K, i, j);
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
			rowvalue = matrix_get(model->V, m+1, K-1, i, j);
			value += pow(rowvalue, 2.0);
		}
	}
	loss += model->lambda * value;

	return loss;
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
			matrix_mul(model->V, m+1, K-1, i, j, 2.0);
			value = - matrix_get(model->Vbar, m+1, K-1, i, j);
			matrix_add(model->V, m+1, K-1, i, j, value);
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
			q = matrix_get(model->Q, model->n, model->K, i, j);
			value = 0.0;
			if (q <= -model->kappa) {
				value = 1.0 - q - (model->kappa+1.0)/2.0;
			} else if (q <= 1.0) {
				value = 1.0/(2.0*model->kappa+2.0)*pow(1.0 - q,
					       	2.0);
			}
			matrix_set(model->H, model->n, model->K, i, j, value);
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
	double q, *uu_row = NULL;

	long n = model->n;
	long K = model->K;

	gensvm_calculate_ZV(model, data, ZV);

	for (i=0; i<n; i++) {
		for (j=0; j<K; j++) {
			if (j == (data->y[i]-1))
				continue;
			uu_row = &matrix_get(model->UU, K*K, K-1, 
					(data->y[i]-1)*K+j, 0);
			#if MAJOR_ORDER == 'r'
			q = cblas_ddot(K-1, &ZV[i*(K-1)], 1, uu_row, 1);
			#else
			q = cblas_ddot(K-1, &ZV[i], n, uu_row, K*K);
			#endif
			matrix_set(model->Q, n, K, i, j, q);
		}
	}
}
