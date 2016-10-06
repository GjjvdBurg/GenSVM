/**
 * @file gensvm_task.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Functions for initializing and freeing a GenTask
 *
 */

#include "gensvm_task.h"

/**
 * @brief Initialize a GenTask structure
 *
 * @details
 * A GenTask structure is initialized and the default value for the
 * parameters are set. A pointer to the initialized GenTask is returned.
 *
 * @returns 	initialized GenTask
 */
struct GenTask *gensvm_init_task()
{
	struct GenTask *t = Malloc(struct GenTask, 1);

	t->kerneltype = K_LINEAR;
	t->weight_idx = 1;
	t->folds = 10;
	t->ID = -1;
	t->p = 1.0;
	t->kappa = 0.0;
	t->lambda = 1.0;
	t->epsilon = 1e-6;
	t->kernelparam = NULL;
	t->train_data = NULL;
	t->test_data = NULL;
	t->performance = 0.0;

	return t;
}

/**
 * @brief Free the GenTask struct
 *
 * @details
 * Freeing the allocated memory of the GenTask means freeing _only_ the 
 * kernelparam array, and the task itself. The datasets are not freed, as 
 * these are shared between all tasks.
 *
 * @param[in] 	t 	GenTask to be freed
 *
 */
void gensvm_free_task(struct GenTask *t)
{
	free(t->kernelparam);
	free(t);
	t = NULL;
}
