/**
 * @file msvmmaj.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Definitions for common structures
 *
 * @details
 * Contains documentation and declarations of MajModel and MajData.
 *
 */

#ifndef MSVMMAJ_H
#define MSVMMAJ_H

#include "globals.h"
#include "types.h"

/**
 * @brief A structure to represent a single MSVMMaj model.
 *
 * @param weight_idx 		which weights to use (1 = unit, 2 = group)
 * @param K 			number of classes in the dataset
 * @param n  			number of instances in the dataset
 * @param m 		 	number of predictors in the dataset
 * @param epsilon 		stopping criterion
 * @param p 			parameter for the L_p norm
 * @param kappa 		parameter for the Huber hinge
 * @param lambda 		regularization parameter
 * @param *W 			pointer to the weight matrix
 * @param *t 			pointer to the translation vector
 * @param *V 	 		pointer to the augmented weight matrix
 * @param *Vbar 		pointer to the augmented weight matrix from a 
 * 				previous iteration
 * @param *U 			pointer to the simplex matrix
 * @param *UU 			pointer to the 3D simplex difference matrix
 * @param *Q 			pointer to the error matrix
 * @param *H 			pointer to the Huber weighted error matrix
 * @param *R 			pointer to the 0-1 auxiliary matrix
 * @param *J 			pointer to the diagonal matrix in the 
 * regularization term
 * @param *rho 			pointer to the instance weight vector
 * @param training_error 	error after training has completed
 * @param *data_file 		pointer to the filename of the data
 * @param kerneltype 		kernel to be used in the model
 * @param kernelparam 		pointer to the vector of kernel parameters
 *
 * @TODO
 * change R to int, it's a binary matrix
 */
struct MajModel {
	int weight_idx; 	
	long K; 	
	long n; 
	long m;
	double epsilon;		
	double p; 	
	double kappa;
	double lambda; 		
	double *W;	
	double *t; 
	double *V; 	
	double *Vbar; 		
	double *U;		
	double *UU; 	
	double *Q;
	double *H; 
	double *R;
	double *J;
	double *rho; 	
	double training_error; 
	char *data_file; 
	KernelType kerneltype; 
	double *kernelparam; 	
};

/**
 * @brief A structure to represent the data.
 *
 * @param K 		number of classes
 * @param n 		number of instances
 * @param m 		number of predictors
 * @param *y 		pointer to vector of class labels
 * @param *Z 		pointer to augmented data matrix
 * @param kerneltype 	kerneltype used in MajData::Z
 * @param *kernelparam 	kernel parameters used in MajData::Z
 *
 */
struct MajData {
	long K; 		
	long n; 		
	long m;		 
 	long *y; 	
	double *Z; 
	KernelType kerneltype; 	
 	double *kernelparam; 
};

#endif
