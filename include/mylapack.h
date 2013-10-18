#ifndef MYLAPACK_H
#define MYLAPACK_H

#include "globals.h"

int dposv(char UPLO, int N, int NRHS, double *A, int LDA, double *B,
		int LDB);
int dsysv(char UPLO, int N, int NRHS, double *A, int LDA, int *IPIV,
		double *B, int LDB, double *WORK, int LWORK);

#endif
