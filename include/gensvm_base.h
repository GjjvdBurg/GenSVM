/**
 * @file gensvm_base.h
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Header file for gensvm_base.c
 *
 * @details
 * Contains documentation and declarations of GenModel and GenData, and 
 * function declarations for the functions in gensvm_base.c.
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

#ifndef GENSVM_BASE_H
#define GENSVM_BASE_H

// includes
#include "gensvm_sparse.h"

// type declarations

/**
 * @brief A structure to represent the data.
 *
 * @param K 		number of classes
 * @param n 		number of instances
 * @param m 		number of predictors
 * @param y 		pointer to vector of class labels
 * @param Z 		pointer to augmented data matrix
 * @param spZ 		pointer to the sparse augmented data matrix
 * @param RAW 		pointer to augmented raw data matrix
 * @param J 		pointer to regularization vector
 * @param Sigma 	eigenvalues from the reduced eigendecomposition
 * @param kerneltype 	kerneltype used in GenData::Z
 * @param gamma 	kernel parameter for RBF, poly, and sigmoid
 * @param coef 		kernel parameter for poly and sigmoid
 * @param degree 	kernel parameter for poly
 *
 */
struct GenData {
	long K;
	///< number of classes
	long n;
	///< number of instances
	long m;
	///< number of predictors (width of RAW)
	long r;
	///< number of eigenvalues (width of Z)
 	long *y;
	///< array of class labels, 1..K
	double *Z;
	///< augmented data matrix (either equal to RAW or to the eigenvectors
	///< of the kernel matrix)
	struct GenSparse *spZ;
	///< sparse representation of the augmented data matrix
	double *RAW;
	///< augmented raw data matrix
	double *Sigma;
	///< eigenvalues from the reduced eigendecomposition
	KernelType kerneltype;
	///< kerneltype used to generate the kernel corresponding to the data 
	///< in Z
	double gamma;
	///< kernel parameter for RBF, poly, and sigmoid
	double coef;
	///< kernel parameter for poly and sigmoid
	double degree;
	///< kernel parameter for poly
};

/**
 * @brief A structure to represent a single GenSVM model.
 *
 */
struct GenModel {
	int weight_idx;
	///< which weights to use (0 = raw, 1 = unit, 2 = group)
	long K;
	///< number of classes in the dataset
	long n;
	///< number of instances in the dataset
	long m;
	///< number of predictor variables in the dataset
	double epsilon;
	///< stopping criterion for the IM algorithm.
	double p;
	///< parameter for the L-p norm in the loss function
	double kappa;
	///< parameter for the Huber hinge function
	double lambda;
	///< regularization parameter in the loss function
	double gamma;
	///< kernel parameter for RBF, poly, and sigmoid
	double coef;
	///< kernel parameter for poly and sigmoid
	double degree;
	///< kernel parameter for poly
	double *V;
	///< augmented weight matrix
	double *Vbar;
	///< augmented weight matrix from the previous iteration of the IM
	///< algorithm
	double *U;
	///< simplex matrix (K x (K-1))
	double *UU;
	///< simplex difference matrix (K*K x (K-1))
	double *Q;
	///< error matrix
	double *H;
	///< Huber weighted error matrix
	double *rho;
	///< vector of instance weights
	double training_error;
	///< loss function value after training has finished
	long elapsed_iter;
	///< number of elapsed iterations in training
	double elapsed_time;
	///< time in seconds elapsed for optimization
	char *data_file;
	///< filename of the data
	KernelType kerneltype;
	///< type of kernel used in the model
	double kernel_eigen_cutoff;
	///< cutoff value for the ratio of eigenvalues in the reduced 
	//eigendecomposition.
	long max_iter;
	///< maximum number of iterations of the algorithm
	int status;
	///< status of the model after training
	long seed;
	///< seed for the random number generator (-1 = random)
};

/**
 * @brief A structure to hold the GenSVM workspace
 *
 */
struct GenWork {
	long n;
	///< number of instances for the workspace
	long m;
	///< number of features for the workspace
	long K;
	///< number of classes for the workspace

	double *LZ;
	///< n x (m+1) working matrix for the Z'*A*Z calculation
	double *ZB;
	///< (m+1) x (K-1) working matrix for the Z'*B calculation
	double *ZBc;
	///< (K-1) x (m+1) working matrix for the Z'*B calculation
	double *ZAZ;
	///< (m+1) x (m+1) working matrix for the Z'*A*Z calculation
	double *tmpZAZ;
	///< (m+1) x (m+1) temporary working matrix for the Z'*A*Z calculation
	double *ZV;
	///< n x (K-1) working matrix for the Z * V calculation
	double *beta;
	///< K-1 working vector for a row of the B matrix

	double *A;
	///< n working vector for the alpha_i values. For CSC sparse only
	double *B;
	///< n x (K-1) working matrix for the beta_i vectors. For CSC sparse
	///< only
};

// function declarations
struct GenModel *gensvm_init_model(void);
void gensvm_allocate_model(struct GenModel *model);
void gensvm_reallocate_model(struct GenModel *model, long n, long m);
void gensvm_free_model(struct GenModel *model);

struct GenData *gensvm_init_data(void);
void gensvm_free_data(struct GenData *data);

struct GenWork *gensvm_init_work(struct GenModel *model);
void gensvm_free_work(struct GenWork *work);
void gensvm_reset_work(struct GenWork *work);

#endif
