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

#ifndef MSVMMAJ_TIMER_H
#define MSVMMAJ_TIMER_H

#include "globals.h"

double elapsed_time(clock_t s_time, clock_t e_time);

#endif
