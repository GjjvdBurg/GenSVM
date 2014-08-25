/**
 * @file gensvm_kernel.h
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Header file for kernel functionality
 *
 * @details
 * Contains function declarations for computing the kernel matrix
 * in nonlinear MSVMGen. Additional kernel functions should be 
 * included here and in gensvm_kernel.c
 *
 */

#ifndef GENSVM_KERNEL_H
#define GENSVM_KERNEL_H

#include "globals.h"

// forward declarations
struct GenData;
struct GenModel;

// function declarations
void gensvm_make_kernel(struct GenModel *model, struct GenData *data);

long gensvm_make_eigen(double *K, long n, double **P, double **Lambda);

void gensvm_make_crosskernel(struct GenModel *model,
	       	struct GenData *data_train, struct GenData *data_test,
		double **K2);

double gensvm_compute_rbf(double *x1, double *x2, double *kernelparam,
		long n);
double gensvm_compute_poly(double *x1, double *x2, double *kernelparam,
		long n);
double gensvm_compute_sigmoid(double *x1, double *x2, double *kernelparam,
		long n);
#endif
