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

// function declarations
double gensvm_elapsed_time(clock_t s_time, clock_t e_time);

#endif
