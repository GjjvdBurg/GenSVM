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

// forward declarations
struct GenData;
struct GenModel;

// function declarations

void gensvm_kernel_preprocess(struct GenModel *model, struct GenData *data);
void gensvm_kernel_postprocess(struct GenModel *model,
	       	struct GenData *traindata, struct GenData *testdata);
void gensvm_make_kernel(struct GenModel *model, struct GenData *data,
		double *K);
long gensvm_make_eigen(double *K, long n, double **P, double **Sigma);
void gensvm_make_crosskernel(struct GenModel *model,
		struct GenData *data_train, struct GenData *data_test,
		double **K2);
void gensvm_make_trainfactor(struct GenData *data, double *P, double *Sigma,
		long r);
void gensvm_make_testfactor(struct GenData *testdata,
	       	struct GenData *traindata, double *K2);
double gensvm_dot_rbf(double *x1, double *x2, double *kernelparam, long n);
double gensvm_dot_poly(double *x1, double *x2, double *kernelparam, long n);
double gensvm_dot_sigmoid(double *x1, double *x2, double *kernelparam, long n);

#endif
