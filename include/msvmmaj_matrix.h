/**
 * @file msvmmaj_matrix.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Header file for msvmmaj_matrix.c
 *
 * @details
 * Contains function declarations for functions useful for dealing with matrices.
 *
 */

#ifndef MSVMMAJ_MATRIX_H
#define MSVMMAJ_MATRIX_H

#include "globals.h"

void matrix_set(double *M, long cols, long i, long j, double val);
void matrix_add(double *M, long cols, long i, long j, double val);
void matrix_mul(double *M, long cols, long i, long j, double val);

double matrix_get(double *M, long cols, long i, long j);

void matrix3_set(double *M, long N2, long N3, long i, long j, long k,
		double val);
double matrix3_get(double *M, long N2, long N3, long i, long j, long k);

void print_matrix(double *M, long rows, long cols);

#endif
