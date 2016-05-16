/**
 * @file gensvm_copy.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Function for copying a GenModel instance
 *
 */

#include "gensvm_copy.h"

/**
 * @brief Copy model parameters between two GenModel structs
 *
 * @details
 * The parameters copied are GenModel::weight_idx, GenModel::epsilon,
 * GenModel::p, GenModel::kappa, and GenModel::lambda.
 *
 * @param[in] 		from 	GenModel to copy parameters from
 * @param[in,out] 	to 	GenModel to copy parameters to
 */
void gensvm_copy_model(struct GenModel *from, struct GenModel *to)
{
	to->weight_idx = from->weight_idx;
	to->epsilon = from->epsilon;
	to->p = from->p;
	to->kappa = from->kappa;
	to->lambda = from->lambda;

	to->kerneltype = from->kerneltype;
	switch (to->kerneltype) {
		case K_LINEAR:
			break;
		case K_POLY:
			to->kernelparam = Malloc(double, 3);
			to->kernelparam[0] = from->kernelparam[0];
			to->kernelparam[1] = from->kernelparam[1];
			to->kernelparam[2] = from->kernelparam[2];
			break;
		case K_RBF:
			to->kernelparam = Malloc(double, 1);
			to->kernelparam[0] = from->kernelparam[0];
			break;
		case K_SIGMOID:
			to->kernelparam = Malloc(double, 2);
			to->kernelparam[0] = from->kernelparam[0];
			to->kernelparam[1] = from->kernelparam[1];
			break;
	}
}
