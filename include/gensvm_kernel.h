/**
 * @file gensvm_kernel.h
 * @author G.J.J. van den Burg
 * @date 2014-01-07
 * @brief Header file for gensvm_kernel.c
 *
 * @details
 * Contains function declarations for computing the kernel matrix
 * in nonlinear MSVMGen. Additional kernel functions should be
 * included here and in gensvm_kernel.c
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

#ifndef GENSVM_KERNEL_H
#define GENSVM_KERNEL_H

// includes
#include "gensvm_base.h"

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
int dsyevx(char JOBZ, char RANGE, char UPLO, int N, double *A, int LDA,
	       	double VL, double VU, int IL, int IU, double ABSTOL,
		int *M, double *W, double *Z, int LDZ, double *WORK, int LWORK,
		int *IWORK, int *IFAIL);
double dlamch(char CMACH);
#endif
