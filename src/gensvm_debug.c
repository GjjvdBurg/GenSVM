/**
 * @file gensvm_debug.c
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Functions facilitating debugging
 *
 * @details
 * Defines functions useful for debugging matrices.
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

#include "gensvm_debug.h"

/**
 * @brief Print a dense matrix
 *
 * @details
 * Debug function to print a matrix
 *
 * @param[in] 	M 	matrix
 * @param[in] 	rows 	number of rows of M
 * @param[in] 	cols 	number of columns of M
 */
void gensvm_print_matrix(double *M, long rows, long cols)
{
	long i, j;

	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			if (j > 0)
				note(" ");
			note("%+6.6f", matrix_get(M, rows, cols, i, j));
		}
		note("\n");
	}
	note("\n");
}

/**
 * @brief Print a sparse matrix
 *
 * @details
 * Debug function to print a GenSparse sparse matrix
 *
 * @param[in] 	A 	a GenSparse matrix to print
 *
 */
void gensvm_print_sparse(struct GenSparse *A)
{
	long i, ix_len;

	// print matrix dimensions
	note("Sparse Matrix:\n");
	note("\ttype = %s\n", ((A->type == CSR) ? "CSR" : "CSC"));
	note("\tnnz = %li, rows = %li, cols = %li\n", A->nnz, A->n_row,
			A->n_col);

	// print nonzero values
	note("\tvalues = [ ");
	for (i=0; i<A->nnz; i++) {
		if (i != 0) note(", ");
		note("%f", A->values[i]);
	}
	note(" ]\n");

	// print cumulative lengths
	note("\tIX = [ ");
	ix_len = (A->type == CSR) ? A->n_row + 1 : A->n_col + 1;
	for (i=0; i<ix_len; i++) {
		if (i != 0) note(", ");
		note("%i", A->ix[i]);
	}
	note(" ]\n");

	// print indices
	note("\tJX = [ ");
	for (i=0; i<A->nnz; i++) {
		if (i != 0) note(", ");
		note("%i", A->jx[i]);
	}
	note(" ]\n");
}


/**
 * @brief Print a GenData structure
 *
 * @param[in] 	data 	GenData structure to print
 *
 */
void gensvm_print_data(struct GenData *data)
{
	char kernel_names[4][8] = {"linear", "poly", "rbf", "sigmoid"};

	note("GenData structure\n");
	note("-----------------\n");
	note("Address: %p\n", data);
	note("\n");
	note("n = %li\n", data->n);
	note("m = %li\n", data->m);
	note("K = %li\n", data->K);
	note("r = %li\n", data->r);
	note("Kernel parameters:\n");
	note("\ttype = %s\n", kernel_names[data->kerneltype]);
	note("\tgamma = %.16f\n", data->gamma);
	note("\tcoef = %.16f\n", data->coef);
	note("\tdegree = %.16f\n", data->degree);
	note("\n");

	note("y:\n");
	if (data->y != NULL) {
		int i;
		for (i=0; i<data->n; i++)
			note("%i ", data->y[i]);
		note("\n");
	}

	if (data->Sigma != NULL) {
		note("Sigma:\n");
		gensvm_print_matrix(data->Sigma, 1, data->r);
	}
	if (data->Z == NULL && data->RAW == NULL) {
		note("spZ:\n");
		gensvm_print_sparse(data->spZ);
	} else {
		note("Z:\n");
		gensvm_print_matrix(data->Z, data->n, data->r+1);
		if (data->Z != data->RAW) {
			note("\nRAW:\n");
			gensvm_print_matrix(data->RAW, data->n, data->m+1);
		}
	}
}

/**
 * @brief Print a GenModel structure
 *
 * @param[in] 	model 	GenModel structure to print
 *
 */
void gensvm_print_model(struct GenModel *model)
{
	char kernel_names[4][8] = {"linear", "poly", "rbf", "sigmoid"};

	note("GenModel structure\n");
	note("------------------\n");
	note("Address: %p\n", model);
	note("Data file: %s\n", model->data_file);
	note("\n");
	note("n = %li\n", model->n);
	note("m = %li\n", model->m);
	note("K = %li\n", model->K);
	note("weight_idx = %i\n", model->weight_idx);
	note("epsilon = %g\n", model->epsilon);
	note("p = %.16f\n", model->p);
	note("kappa = %.16f\n", model->kappa);
	note("lambda = %.16f\n", model->lambda);
	note("max_iter = %li\n", model->max_iter);
	note("seed = %li\n", model->seed);
	note("Kernel parameters:\n");
	note("\ttype = %s\n", kernel_names[model->kerneltype]);
	note("\tgamma = %.16f\n", model->gamma);
	note("\tcoef = %.16f\n", model->coef);
	note("\tdegree = %.16f\n", model->degree);
	note("\tkernel_eigen_cutoff = %.16f\n", model->kernel_eigen_cutoff);
	note("Results:\n");
	note("\ttraining_error = %.16f\n", model->training_error);
	note("\telapsed_iter = %li\n", model->elapsed_iter);
	note("\tstatus = %i\n", model->status);

	note("\nV:\n");
	if (model->V != NULL)
		gensvm_print_matrix(model->V, model->m+1, model->K-1);

	note("\nVbar:\n");
	if (model->Vbar != NULL)
		gensvm_print_matrix(model->Vbar, model->m+1, model->K-1);

	note("\nU:\n");
	if (model->U != NULL)
		gensvm_print_matrix(model->U, model->K, model->K-1);

	note("\nUU:\n");
	if (model->UU != NULL)
		gensvm_print_matrix(model->UU, model->K * model->K, model->K-1);

	note("\nQ:\n");
	if (model->Q != NULL)
		gensvm_print_matrix(model->Q, model->n, model->K);

	note("\nH:\n");
	if (model->H != NULL)
		gensvm_print_matrix(model->H, model->n, model->K);

	note("\nrho:\n");
	if (model->rho != NULL)
		gensvm_print_matrix(model->rho, 1, model->n);
}
