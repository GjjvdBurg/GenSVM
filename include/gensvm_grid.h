/**
 * @file gensvm_grid.h
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Header file for gensvm_grid.c
 *
 * @details
 * The grid search for the optimal parameters is done through a queue.
 * This file contains struct definitions for this queue and a single
 * task in a queue, as well as a structure for the complete training
 * scheme. Function declarations are also included.
 *
 * @copyright
 Copyright 2016, G.J.J. van den Burg.

 This file is part of GenSVM.

 GenSVM is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 GenSVM is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with GenSVM. If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef GENSVM_GRID_H
#define GENSVM_GRID_H

#include "gensvm_globals.h"

/**
 * @brief Structure for describing the entire grid search
 *
 * @param traintype 		type of training to use
 * @param kerneltype 		type of kernel to use throughout training
 * @param repeats 		number of repeats to be done after the grid
 * 				search to find the parameter set with the
 * 				most consistent high performance
 * @param folds 		number of folds in cross validation
 * @param Np 			size of the array of p values
 * @param Nl 			size of the array of lambda values
 * @param Nk 			size of the array of kappa values
 * @param Ne 			size of the array of epsilon values
 * @param Nw 			size of the array of weight_idx values
 * @param Ng 			size of the array of gamma values
 * @param Nc 			size of the array of coef values
 * @param Nd 			size of the array of degree values
 * @param *weight_idxs 		array of weight_idxs
 * @param *ps 			array of p values
 * @param *lambdas 		array of lambda values
 * @param *kappas 		array of kappa values
 * @param *epsilons 		array of epsilon values
 * @param *gammas 		array of gamma values
 * @param *coefs 		array of coef values
 * @param *degrees 		array of degree values
 * @param *train_data_file 	filename of train data file
 * @param *test_data_file 	filename of test data file
 *
 */
struct GenGrid {
	TrainType traintype;
	///< type of training to use
	KernelType kerneltype;
	///< type of kernel to use throughout training
	long folds;
	///< number of folds in cross validation
	long repeats;
	///< number of repeats to be done after the grid search to find the
	///< parameter set with the most consistent high performance
	double percentile;
	///< percentile to use for the consistency repeats
	long Np;
	///< size of the array of p values
	long Nl;
	///< size of the array of lambda values
	long Nk;
	///< size of the array of kappa values
	long Ne;
	///< size of the array of epsilon values
	long Nw;
	///< size of the array of weight_idx values
	long Ng;
	///< size of the array of gamma values
	long Nc;
	///< size of the array of coef values
	long Nd;
	///< size of the array of degree values
	int *weight_idxs;
	///< array of weight_idxs
	double *ps;
	///< array of p values
	double *lambdas;
	///< array of lambda values
	double *kappas;
	///< array of kappa values
	double *epsilons;
	///< array of epsilon values
	double *gammas;
	///< array of gamma values
	double *coefs;
	///< array of coef values
	double *degrees;
	///< array of degree values
	char *train_data_file;
	///< filename of train data file
	char *test_data_file;
	///< filename of test data file
};

// function declarations
struct GenGrid *gensvm_init_grid(void);
void gensvm_free_grid(struct GenGrid *grid);

#endif
