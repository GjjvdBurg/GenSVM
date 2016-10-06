/**
 * @file gensvm_queue.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Functions for initializing and freeing a GenQueue
 *
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
struct GenQueue *gensvm_init_queue()
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
