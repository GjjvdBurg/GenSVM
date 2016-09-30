/**
 * @file gensvm_base.h
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Header file for gensvm_base.c
 *
 * @details
 * Contains documentation and declarations of GenModel and GenData, and 
 * function declarations for the functions in gensvm_base.c.
 *
 */

#ifndef GENSVM_BASE_H
#define GENSVM_BASE_H

// includes
#include "gensvm_globals.h"

// type declarations

/**
 * @brief A structure to represent the data.
 *
 * @param K 		number of classes
 * @param n 		number of instances
 * @param m 		number of predictors
 * @param *y 		pointer to vector of class labels
 * @param *Z 		pointer to augmented data matrix
 * @param *RAW 		pointer to augmented raw data matrix
 * @param *J 		pointer to regularization vector
 * @param kerneltype 	kerneltype used in GenData::Z
 * @param *kernelparam 	kernel parameters used in GenData::Z
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
	double *RAW;
	///< augmented raw data matrix
	double *Sigma;
	KernelType kerneltype;
 	double *kernelparam;
};

/**
 * @brief A structure to represent a single GenSVM model.
 *
 */
struct GenModel {
	int weight_idx;
	///< which weights to use (1 = unit, 2 = group)
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
	double *V;
	///< augmented weight matrix
	double *Vbar;
	///< augmented weight matrix from the previous iteration of the IM
	///< algorithm
	double *U;
	///< simplex matrix
	double *UU;
	///< simplex difference matrix
	double *Q;
	///< error matrix
	double *H;
	///< Huber weighted error matrix
	double *rho;
	///< vector of instance weights
	double training_error;
	///< loss function value after training has finished
	char *data_file;
	///< filename of the data
	KernelType kerneltype;
	///< type of kernel used in the model
	double *kernelparam;
	///< array of kernel parameters, size depends on kernel type
};

// function declarations
struct GenModel *gensvm_init_model();
void gensvm_allocate_model(struct GenModel *model);
void gensvm_reallocate_model(struct GenModel *model, long n, long m);
void gensvm_free_model(struct GenModel *model);
struct GenData *gensvm_init_data();
void gensvm_free_data(struct GenData *data);

#endif
