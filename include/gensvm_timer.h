/**
 * @file timer.h
 * @author Gertjan van den Burg
 * @date August, 2013
 * @brief Header file for timer.c
 *
 * @details
 * Function declaration for timer function used to measure computation time.
 *
 */

#ifndef GENSVM_TIMER_H
#define GENSVM_TIMER_H

// includes
#include "gensvm_globals.h"

/// Timer macro for easily recording time
#define Timer(spec) clock_gettime(CLOCK_MONOTONIC_RAW, &spec)

// function declarations
double gensvm_elapsed_time(struct timespec *start, struct timespec *stop);

#endif
