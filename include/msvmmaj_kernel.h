/**
 * @file msvmmaj_kernel.h
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Header file for kernel functionality
 *
 * @details
 * Contains function declarations for computing the kernel matrix
 * in nonlinear MSVMMaj. Additional kernel functions should be 
 * included here and in msvmmaj_kernel.c
 *
 */

#ifndef MSVMMAJ_KERNEL_H
#define MSVMMAJ_KERNEL_H

#include "globals.h"

// forward declarations
struct MajData;
struct MajModel;

// function declarations
void msvmmaj_make_kernel(struct MajModel *model, struct MajData *data);

double msvmmaj_compute_rbf(double *x1, double *x2, double *kernelparam,
		long n);
double msvmmaj_compute_poly(double *x1, double *x2, double *kernelparam,
		long n);
double msvmmaj_compute_sigmoid(double *x1, double *x2, double *kernelparam,
		long n);
#endif
