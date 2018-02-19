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
void gensvm_kernel_copy_kernelparam_to_data(struct GenModel *model, 
		struct GenData *data);
void gensvm_kernel_preprocess(struct GenModel *model, struct GenData *data);
void gensvm_kernel_postprocess(struct GenModel *model,
	       	struct GenData *traindata, struct GenData *testdata);
void gensvm_kernel_compute(struct GenModel *model, struct GenData *data,
		double *K);
long gensvm_kernel_eigendecomp(double *K, long n, double cutoff, 
		double **P_ret, double **Sigma_ret);
double *gensvm_kernel_cross(struct GenModel *model, struct GenData *data_train,
		struct GenData *data_test);
void gensvm_kernel_trainfactor(struct GenData *data, double *P, double *Sigma,
		long r);
void gensvm_kernel_testfactor(struct GenData *testdata,
	       	struct GenData *traindata, double *K2);

double gensvm_kernel_dot_rbf(double *x, double *y, long n, long incx, 
		long incy, double gamma);
double gensvm_kernel_dot_poly(double *x, double *y, long n, long incx, 
		long incy, double gamma, double coef, double degree);
double gensvm_kernel_dot_sigmoid(double *x, double *y, long n, long incx, 
		long incy, double gamma, double coef);

int dsyevx(char JOBZ, char RANGE, char UPLO, int N, double *A, int LDA,
	       	double VL, double VU, int IL, int IU, double ABSTOL,
		int *M, double *W, double *Z, int LDZ, double *WORK, int LWORK,
		int *IWORK, int *IFAIL);
double dlamch(char CMACH);
#endif
