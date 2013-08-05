#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cblas.h>
#include <string.h>
#include "util.h"

void simplex_gen(long K, double *U);
void category_matrix(struct Model *model, struct Data *data);
void simplex_diff(struct Model *model, struct Data *dataset);

void calculate_errors(struct Model *model, struct Data *data, double *ZV);
void calculate_huber(struct Model *model);

double get_msvmmaj_loss(struct Model *model, struct Data *data, double *ZV);
void msvmmaj_update(struct Model *model, struct Data *data, 
		int *ClassIdx, double *A, double *B, double *Omega, 
		double *ZAZ, double *ZAZV, double *ZAZVT);
void step_doubling(struct Model *model);

void main_loop(struct Model *model, struct Data *data);

int dposv(char UPLO, int N, int NRHS, double *A, int LDA, double *B, int LDB);

void initialize_weights(struct Data *data, struct Model *model);
