/**
 * @file gensvm_lapack.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Header file for gensvm_lapack.c
 * 
 * @details
 * Function declarations for external LAPACK functions 
 *
 */

#ifndef GENSVM_LAPACK_H
#define GENSVM_LAPACK_H

#include "globals.h"

int dposv(char UPLO, int N, int NRHS, double *A, int LDA, double *B,
		int LDB);
int dsysv(char UPLO, int N, int NRHS, double *A, int LDA, int *IPIV,
		double *B, int LDB, double *WORK, int LWORK);
int dsyevx(char JOBZ, char RANGE, char UPLO, int N, double *A, int LDA,
	       	double VL, double VU, int IL, int IU, double ABSTOL,
		int *M, double *W, double *Z, int LDZ, double *WORK, int LWORK,
		int *IWORK, int *IFAIL);
double dlamch(char CMACH);
#endif
