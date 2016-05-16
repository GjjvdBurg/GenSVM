/**
 * @file gensvm_queue.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Header file for gensvm_queue.c
 *
 * @details
 * The grid search for the optimal parameters is done through a queue.
 * This file contains struct definitions for this queue. Function declarations 
 * for initializing and freeing the queue are also included.
 *
 */

#ifndef GENSVM_QUEUE_H
#define GENSVM_QUEUE_H

#include "gensvm_task.h"

/**
 * @brief Simple task queue.
 *
 * This struct is basically just an array of pointers to Task instances,
 * with a length and an index of the current task.
 *
 * @param **tasks 	array of pointers to Task structs
 * @param N 		size of task array
 * @param i 		index used for keeping track of the queue
 */
struct GenQueue {
	struct GenTask **tasks;
	long N;
	long i;
};

// function declarations
struct GenQueue *gensvm_init_queue();
void gensvm_free_queue(struct GenQueue *q);
struct GenTask *get_next_task(struct GenQueue *q);

#endif
