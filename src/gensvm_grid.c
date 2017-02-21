/**
 * @file gensvm_grid.c
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Functions for initializing GenGrid structures
 *
 * @details
 * This file contains functions for initializing and freeing a GenGrid 
 * instance.  In addition, default values for this structure are defined here 
 * (and only here).
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

#include "gensvm_grid.h"

/**
 * @brief Initialize a GenGrid structure
 *
 * @brief
 * This function is used to initialize a GenGrid struct, and set its default 
 * parameters. A pointer to the generated struct is returned.
 *
 * @return  	initialized GenGrid struct
 *
 */
struct GenGrid *gensvm_init_grid(void)
{
	struct GenGrid *grid = Malloc(struct GenGrid, 1);

	// initialize to defaults
	grid->traintype = CV;
	grid->kerneltype = K_LINEAR;
	grid->folds = 10;
	grid->repeats = 0;
	grid->percentile = 95.0;
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

/**
 * @brief Free a GenGrid structure
 *
 * @details
 * This function frees all elements of a GenGrid structure, including the 
 * GenGrid structure itself. The provided argument is set to NULL on exit.
 *
 * @param[in] 	grid 	a GenGrid struct to free
 *
 */
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
