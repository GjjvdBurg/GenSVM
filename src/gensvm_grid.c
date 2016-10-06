/**
 * @file gensvm_grid.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Functions for initializing GenGrid structures
 *
 * @details
 * This file contains functions for initializing and freeing a GenGrid 
 * instance.  In addition, default values for this structure are defined here 
 * (and only here).
 *
 */

#include "gensvm_grid.h"

struct GenGrid *gensvm_init_grid()
{
	struct GenGrid *grid = Malloc(struct GenGrid, 1);

	// initialize to defaults
	grid->traintype = CV;
	grid->kerneltype = K_LINEAR;
	grid->repeats = 0;
	grid->folds = 10;
	grid->Np = 0;
	grid->Nl = 0;
	grid->Nk = 0;
	grid->Ne = 0;
	grid->Nw = 0;
	grid->Ng = 0;
	grid->Nc = 0;
	grid->Nd = 0;

	// set arrays to NULL
	grid->weight_idxs = NULL;
	grid->ps = NULL;
	grid->lambdas = NULL;
	grid->kappas = NULL;
	grid->epsilons = NULL;
	grid->gammas = NULL;
	grid->coefs = NULL;
	grid->degrees = NULL;
	grid->train_data_file = NULL;
	grid->test_data_file = NULL;

	return grid;
}

void gensvm_free_grid(struct GenGrid *grid)
{
	free(grid->weight_idxs);
	free(grid->ps);
	free(grid->lambdas);
	free(grid->kappas);
	free(grid->epsilons);
	free(grid->gammas);
	free(grid->coefs);
	free(grid->degrees);
	free(grid->train_data_file);
	free(grid->test_data_file);
	free(grid);
	grid = NULL;
}
