/**
 * @file msvmmaj_lapack.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Header file for msvmmaj_lapack.c
 * 
 * @details
 * Function declarations for external LAPACK functions 
 *
 */

#ifndef MSVMMAJ_LAPACK_H
#define MSVMMAJ_LAPACK_H

#include "globals.h"

int dposv(char UPLO, int N, int NRHS, double *A, int LDA, double *B,
		int LDB);
int dsysv(char UPLO, int N, int NRHS, double *A, int LDA, int *IPIV,
		double *B, int LDB, double *WORK, int LWORK);

#endif
