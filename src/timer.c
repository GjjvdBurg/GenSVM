/**
 * @file timer.c
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Function for calculating time difference
 *
 * @details
 * This file contains a simple function for calculating the time in seconds
 * elapsed between two clock() calls.
 */

#include <time.h>

#include "timer.h"

/**
 * @brief Calculate the time between two clocks
 *
 * @param[in] 	s_time 	starting time
 * @param[in] 	e_time 	end time
 * @returns 		time elapsed in seconds
 */
double elapsed_time(clock_t s_time, clock_t e_time)
{
	return ((double) (e_time - s_time))/((double) CLOCKS_PER_SEC);
}
