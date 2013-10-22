#ifndef MSVMMAJ_H
#define MSVMMAJ_H

#include "globals.h"
#include "types.h"

/*
	Model structure
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
	double *rho;
	double training_error;
	char *data_file;
	KernelType kerneltype;
	double *kernelparam;
};

/*
	Data structure
*/
struct MajData {
	long K;
	long n;
	long m;
	long *y;
	double *Z;
};

#endif
