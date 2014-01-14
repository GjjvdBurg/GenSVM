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
#include "msvmmaj_matrix.h"
#include "msvmmaj_pred.h"

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
void msvmmaj_predict_labels(struct MajData *data, struct MajModel *model, long *predy)
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
				S[k] = matrix_get(ZV, K-1, i, k) - matrix_get(U, K-1, j, k);
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
