/**
 * @file mylapack.c
 * @author Gertjan van den Burg (burg@ese.eur.nl)
 * @date August 9, 2013
 * @brief Utility functions for interacting with LAPACK
 *
 * @details
 * Functions in this file are auxiliary functions which make it easier
 * to use LAPACK functions from liblapack.
 */

#include "mylapack.h"

/**
 * @name dposv
 * @brief Solve a system of equations AX = B where A is symmetric positive definite.
 * @ingroup libMSVMMaj
 *
 * See the LAPACK documentation at: 
 * http://www.netlib.org/lapack/explore-html/dc/de9/group__double_p_osolve.html
 */
int dposv(char UPLO, int N, int NRHS, double *A, int LDA, double *B,
		int LDB)
{
	extern void dposv_(char *UPLO, int *Np, int *NRHSp, double *A,
			int *LDAp, double *B, int *LDBp, int *INFOp);
	int INFO;
	dposv_(&UPLO, &N, &NRHS, A, &LDA, B, &LDB, &INFO);
	return INFO;
}

/**
 * @name dsysv
 * @brief Solve a system of equations AX = B where A is symmetric.
 * @ingroup libMSVMMaj
 *
 * See the LAPACK documentation at:
 * http://www.netlib.org/lapack/explore-html/d6/d0e/group__double_s_ysolve.html
 */
int dsysv(char UPLO, int N, int NRHS, double *A, int LDA, int *IPIV,
		double *B, int LDB, double *WORK, int LWORK)
{
	extern void dsysv_(char *UPLO, int *Np, int *NRHSp, double *A,
			int *LDAp, int *IPIV, double *B, int *LDBp,
			double *WORK, int *LWORK, int *INFOp);
	int INFO;
	dsysv_(&UPLO, &N, &NRHS, A, &LDA, IPIV, B, &LDB, WORK, &LWORK, &INFO);
	return INFO;
}
