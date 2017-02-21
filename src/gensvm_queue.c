/**
 * @file gensvm_queue.c
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Functions for initializing and freeing a GenQueue
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

#include "gensvm_queue.h"

/**
 * @brief Initialize a GenQueue structure
 *
 * @details
 * A GenQueue structure is initialized and the default value for the
 * parameters are set. A pointer to the initialized queue  is returned.
 *
 * @returns 	initialized GenQueue
 */
struct GenQueue *gensvm_init_queue(void)
{
	struct GenQueue *q = Malloc(struct GenQueue, 1);

	q->tasks = NULL;
	q->N = 0;
	q->i = 0;

	return q;
}

/**
 * @brief Free the GenQueue struct
 *
 * @details
 * Freeing the allocated memory of the GenQueue means freeing every GenTask
 * struct and then freeing the Queue.
 *
 * @param[in] 	q 	GenQueue to be freed
 *
 */
void gensvm_free_queue(struct GenQueue *q)
{
	long i;
	for (i=0; i<q->N; i++) {
		gensvm_free_task(q->tasks[i]);
	}
	free(q->tasks);
	free(q);
	q = NULL;
}

/**
 * @brief Get new GenTask from GenQueue
 *
 * @details
 * Return a pointer to the next GenTask in the GenQueue. If no GenTask 
 * instances are left, NULL is returned. The internal counter GenQueue::i is 
 * used for finding the next GenTask.
 *
 * @param[in] 	q 	GenQueue instance
 * @returns 		pointer to next GenTask
 *
 */
struct GenTask *get_next_task(struct GenQueue *q)
{
	long i = q->i;
	if (i < q->N) {
		q->i++;
		return q->tasks[i];
	}
	return NULL;
}
