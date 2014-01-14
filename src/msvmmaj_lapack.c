/**
 * @file msvmmaj_lapack.c
 * @author Gertjan van den Burg
 * @date August 9, 2013
 * @brief Utility functions for interacting with LAPACK
 *
 * @details
 * Functions in this file are auxiliary functions which make it easier
 * to use LAPACK functions from liblapack.
 */

#include "msvmmaj_lapack.h"

/**
 * @brief Solve AX = B where A is symmetric positive definite.
 *
 * @details
 * Solve a linear system of equations AX = B where A is symmetric positive 
 * definite. This function uses the externel LAPACK routine dposv.
 *
 * @param[in] 		UPLO 	which triangle of A is stored
 * @param[in] 		N 	order of A
 * @param[in] 		NRHS 	number of columns of B
 * @param[in,out] 	A 	double precision array of size (LDA, N). On 
 * 				exit contains the upper or lower factor of the
 * 				Cholesky factorization of A.
 * @param[in] 		LDA 	leading dimension of A
 * @param[in,out] 	B 	double precision array of size (LDB, NRHS). On
 * 				exit contains the N-by-NRHS solution matrix X.
 * @param[in] 		LDB 	the leading dimension of B
 * @returns 			info parameter which contains the status of the
 * 				computation:
 * 					- =0: 	success
 * 					- <0: 	if -i, the i-th argument had
 * 						an illegal value
 * 					- >0: 	if i, the leading minor of A
 * 						was not positive definite
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
 * @brief Solve a system of equations AX = B where A is symmetric.
 *
 * @details
 * Solve a linear system of equations AX = B where A is symmetric. This 
 * function uses the external LAPACK routine dsysv.
 *
 * @param[in] 		UPLO 	which triangle of A is stored
 * @param[in] 		N 	order of A
 * @param[in] 		NRHS 	number of columns of B
 * @param[in,out] 	A 	double precision array of size (LDA, N). On 
 * 				exit contains the block diagonal matrix D and
 * 				the multipliers used to obtain the factor U or
 * 				L from the factorization A = U*D*U**T or
 * 				A = L*D*L**T.
 * @param[in] 		LDA 	leading dimension of A
 * @param[in] 		IPIV 	integer array containing the details of D
 * @param[in,out] 	B 	double precision array of size (LDB, NRHS). On
 * 				exit contains the N-by-NRHS matrix X
 * @param[in] 		LDB 	leading dimension of B
 * @param[out] 		WORK 	double precision array of size max(1,LWORK). On
 * 				exit, WORK(1) contains the optimal LWORK
 * @param[in] 		LWORK 	the length of WORK, can be used for determining
 * 				the optimal blocksize for dsystrf.
 * @returns 			info parameter which contains the status of the
 * 				computation:
 * 					- =0: 	success
 * 					- <0: 	if -i, the i-th argument had an
 * 						illegal value
 * 					- >0: 	if i, D(i, i) is exactly zero,
 * 						no solution can be computed.
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

/**
 * @brief Compute the Cholesky factorization of a real symmetric positive 
 *  	  definite matrix.
 *
 * @details
 * This function uses the external LAPACK routine dpotrf.
 *
 * @param[in] 		UPLO 	which triangle of A is stored
 * @param[in] 		N 	order of A
 * @param[in,out] 	A 	double precision array of size (LDA, N). On
 * 				exit contains the factor U or L of the Cholesky
 * 				factorization
 * @param[in] 		LDA 	leading dimension of A
 * @returns 			info parameter which contains the status of the
 * 				computation:
 * 					- =0: 	success
 * 					- <0: 	if -i, the i-th argument had an
 * 						illegal value
 * 					- >0: 	if i, the leading minor of
 * 						order i is not positive
 * 						definite
 *
 * See the LAPACK documentation at:
 * http://www.netlib.org/lapack/explore-html/d0/d8a/dpotrf_8f.html
 */
int dpotrf(char UPLO, int N, double *A, int LDA)
{
	extern void dpotrf_(char *UPLO, int *N, double *A, int *LDA, int *INFOp);
	int INFO;
	dpotrf_(&UPLO, &N, A, &LDA, &INFO);
	return INFO;
}
