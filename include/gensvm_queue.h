/**
 * @file gensvm_queue.h
 * @author G.J.J. van den Burg
 * @date 2013-08-01
 * @brief Header file for gensvm_queue.c
 *
 * @details
 * The grid search for the optimal parameters is done through a queue.
 * This file contains struct definitions for this queue. Function declarations 
 * for initializing and freeing the queue are also included.
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

#ifndef GENSVM_QUEUE_H
#define GENSVM_QUEUE_H

#include "gensvm_task.h"

/**
 * @brief Simple task queue.
 *
 * This struct is basically just an array of pointers to Task instances,
 * with a length and an index of the current task.
 *
 * @param tasks 	array of pointers to Task structs
 * @param N 		size of task array
 * @param i 		index used for keeping track of the queue
 */
struct GenQueue {
	struct GenTask **tasks;
	///< array of pointers to Task structs
	long N;
	///< size of task array
	long i;
	///< index used for keeping track of the queue
};

// function declarations
struct GenQueue *gensvm_init_queue(void);
void gensvm_free_queue(struct GenQueue *q);
struct GenTask *get_next_task(struct GenQueue *q);

#endif
