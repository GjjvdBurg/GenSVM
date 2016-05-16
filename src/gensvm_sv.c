/**
 * @file gensvm_sv.c
 * @author Gertjan van den Burg
 * @date May, 2014
 * @brief Calculate the number of support vectors
 *
 * @details
 * The function in this file can be used to calculate the number of support
 * vectors are left in a model.
 *
 */

#include "gensvm_sv.h"

/**
 * @brief Calculate the number of support vectors in a model
 *
 * @details
 * If an object is correctly classified, the number of classes for which the
 * error q is larger than 1, is K-1 (i.e., there is no error w.r.t. any of the
 * other classes).  All objects for which this is not the case are thus support
 * vectors.
 *
 * @param[in] 	model 	GenModel with solution
 * @param[in] 	data 	GenData to be used
 * @return 		number of support vectors with this solution
 *
 */
long gensvm_num_sv(struct GenModel *model, struct GenData *data)
{
	long i, j, num_correct, num_sv = 0;
	double value;

	for (i=0; i<data->n; i++) {
		num_correct = 0;
		for (j=0; j<data->K; j++) {
			value = matrix_get(model->Q, data->K, i, j);
			num_correct += (value > 1);
		}
		num_sv += (num_correct < data->K - 1);
	}

	return num_sv;
}
