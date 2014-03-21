/**
 * @file msvmmaj_pred.c
 * @author Gertjan van den Burg
 * @date August 9, 2013
 * @brief Main functions for predicting class labels..
 *
 * @details
 * This file contains functions for predicting the class labels of instances
 * and a function for calculating the predictive performance (hitrate) of 
 * a prediction given true class labels.
 *
 */

#include <cblas.h>

#include "libMSVMMaj.h"
#include "msvmmaj.h"
#include "msvmmaj_kernel.h"
#include "msvmmaj_matrix.h"
#include "msvmmaj_pred.h"

#include "util.h" // testing

void msvmmaj_predict_labels(struct MajData *data_test,
	       	struct MajData *data_train, struct MajModel *model,
	       	long *predy)
{
	if (model->kerneltype == K_LINEAR)
		msvmmaj_predict_labels_linear(data_test, model, predy);
	else
		msvmmaj_predict_labels_kernel(data_test, data_train, model,
			       	predy);
}

/**
 * @brief Predict class labels of data given and output in predy
 *
 * @details
 * The labels are predicted by mapping each instance in data to the 
 * simplex space using the matrix V in the given model. Next, for each
 * instance the nearest simplex vertex is determined using an Euclidean
 * norm. The nearest simplex vertex determines the predicted class label,
 * which is recorded in predy.
 *
 * @param[in] 	data 		MajData to predict labels for
 * @param[in] 	model 		MajModel with optimized V
 * @param[out] 	predy 		pre-allocated vector to record predictions in
 */
void msvmmaj_predict_labels_linear(struct MajData *data,
	       	struct MajModel *model, long *predy)
{
	long i, j, k, label;
	double norm, min_dist;

	long n = data->n; // note that model->n is the size of the training sample.
	long m = data->m;
	long K = model->K; //data->K does not necessarily equal the original K.

	double *S = Calloc(double, K-1);
	double *ZV = Calloc(double, n*(K-1));
	double *U = Calloc(double, K*(K-1));

	// Get the simplex matrix
	msvmmaj_simplex_gen(K, U);

	// Generate the simplex-space vectors
	cblas_dgemm(
			CblasRowMajor,
			CblasNoTrans,
			CblasNoTrans,
			n,
			K-1,
			m+1,
			1.0,
			data->Z,
			m+1,
			model->V,
			K-1,
			0.0,
			ZV,
			K-1);

	// Calculate the distance to each of the vertices of the simplex.
	// The closest vertex defines the class label.
	for (i=0; i<n; i++) {
		label = 0;
		min_dist = 1000000000.0;
		for (j=0; j<K; j++) {
			for (k=0; k<K-1; k++) {
				S[k] = matrix_get(ZV, K-1, i, k) -
				       	matrix_get(U, K-1, j, k);
			}
			norm = cblas_dnrm2(K, S, 1);
			if (norm < min_dist) {
				label = j+1; // labels start counting from 1
				min_dist = norm;
			}
		}
		predy[i] = label;
	}

	free(ZV);
	free(U);
	free(S);
}

void msvmmaj_predict_labels_kernel(struct MajData *data_test,
	       	struct MajData *data_train, struct MajModel *model,
	       	long *predy)
{
	long i, j, k, label;
	double norm, min_dist;

	long n_train = data_train->n;
	long n_test = data_test->n;
	long r = model->m;
	long K = model->K;

	double *K2 = NULL;
	msvmmaj_make_crosskernel(model, data_train, data_test, &K2);

	//printf("K2:\n");
	//print_matrix(K2, n_test, n_train);

	//printf("P:\n");
	//print_matrix(data_train->Z, n_train, r+1);

	//printf("Sigma:\n");
	//print_matrix(data_train->J, 1, r+1);

	double *S = Calloc(double, K-1);
	double *ZV = Calloc(double, n_test*(r+1));
	double *KPS = Calloc(double, n_test*(r+1));
	double *U = Calloc(double, K*(K-1));

	msvmmaj_simplex_gen(K, U);
	
	// were doing the computations explicitly since P is included in 
	// data_train->Z. Might want to look at this some more if it turns out 
	// to be slow.

	double value, rowvalue;
	for (i=0; i<n_test; i++) {
		for (j=1; j<r+1; j++) {
			value = 0.0;
			for (k=0; k<n_train; k++) {
				rowvalue = matrix_get(K2, n_train, i, k);
				rowvalue *= matrix_get(data_train->Z, r+1, k,
					       	j);
				value += rowvalue;
			}
			value *= matrix_get(data_train->J, 1, j, 0);
			matrix_set(KPS, r+1, i, j, value);
		}
		matrix_set(KPS, r+1, i, 0, 1.0);
	}
	
	//printf("\nPrinting KPS:\n");
	//print_matrix(KPS, n_test, r+1);

	//printf("\nPrinting Omega (model->m = %li):\n", model->m);
	//print_matrix(model->V, r+1, K-1);

	cblas_dgemm(
			CblasRowMajor,
			CblasNoTrans,
			CblasNoTrans,
			n_test,
			K-1,
			r+1,
			1.0,
			KPS,
			r+1,
			model->V,
			K-1,
			0.0,
			ZV,
			K-1);

	//printf("\nPrinting ZV:\n");
	//print_matrix(ZV, n_test, K-1);
	//printf("\n");

	for (i=0; i<n_test; i++) {
		label = 0;
		min_dist = 1e10;
		for (j=0; j<K; j++) {
			for (k=0; k<K-1; k++) {
				S[k] = matrix_get(ZV, K-1, i, k) -
					matrix_get(U, K-1, j, k);
			}
			norm = cblas_dnrm2(K, S, 1);
			if (norm < min_dist) {
				label = j+1;
				min_dist = norm;
			}
		}
		predy[i] = label;
		//printf("predy[%li] = %li\n", i, label);
	}

	free(ZV);
	free(U);
	free(S);
	free(KPS);
	free(K2);

	//exit(1);
}

/**
 * @brief Calculate the predictive performance (percentage correct)
 *
 * @details
 * The predictive performance is calculated by simply counting the number
 * of correctly classified samples and dividing by the total number of 
 * samples, multiplying by 100.
 *
 * @param[in] 	data 	the MajData dataset with known labels
 * @param[in] 	predy 	the predicted class labels
 *
 * @returns percentage correctly classified.
 */
double msvmmaj_prediction_perf(struct MajData *data, long *predy)
{
	long i, correct = 0;
	double performance;

	for (i=0; i<data->n; i++)
		if (data->y[i] == predy[i])
			correct++;

	performance = ((double) correct)/((double) data->n)* 100.0;

	return performance;
}
