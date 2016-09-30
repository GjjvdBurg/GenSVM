/**
 * @file gensvm_train.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Header file for gensvm_train.c
 *
 * @details
 * Contains function declarations for functions used to train a single
 * GenModel.
 *
 */

#ifndef GENSVM_OPTIMIZE_H
#define GENSVM_OPTIMIZE_H

#include "gensvm_sv.h"
#include "gensvm_print.h"
#include "gensvm_simplex.h"

// function declarations
void gensvm_optimize(struct GenModel *model, struct GenData *data);
double gensvm_get_loss(struct GenModel *model, struct GenData *data,
		double *ZV);

double gensvm_calculate_omega(struct GenModel *model, long i);
bool gensvm_majorize_is_simple(struct GenModel *model, long i);
void gensvm_calculate_ab_non_simple(struct GenModel *model, long i, long j,
		double *a, double *b_aq);
void gensvm_calculate_ab_simple(struct GenModel *model, long i, long j,
		double *a, double *b_aq);
void gensvm_update_B(struct GenModel *model, long i, long j, double b_aq,
		double omega, double *B);
double gensvm_get_Avalue_update_B(struct GenModel *model, long i, double *B);

void gensvm_get_update(struct GenModel *model, struct GenData *data,
		double *B, double *ZAZ, double *ZAZV, double *ZAZVT);
void gensvm_category_matrix(struct GenModel *model, struct GenData *data);
void gensvm_simplex_diff(struct GenModel *model, struct GenData *dataset);
void gensvm_calculate_errors(struct GenModel *model, struct GenData *data,
		double *ZV);
void gensvm_calculate_huber(struct GenModel *model);
void gensvm_step_doubling(struct GenModel *model);
int dposv(char UPLO, int N, int NRHS, double *A, int LDA, double *B, int LDB);
int dsysv(char UPLO, int N, int NRHS, double *A, int LDA, int *IPIV, double *B,
	       	int LDB, double *WORK, int LWORK);

#endif
