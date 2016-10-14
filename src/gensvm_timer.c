/**
 * @file gensvm_timer.c
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Utility functions relating to time
 *
 * @details
 * This file contains a simple function for calculating the time in seconds
 * elapsed between two Timer() calls.
 *
 */

#include "gensvm_timer.h"

/**
 * @brief Calculate the time between two time recordings
 *
 * @details
 * This function should be used with time recordings done by clock_gettime() 
 * using CLOCK_MONOTONIC_RAW. For this, the Timer() macro has been defined in 
 * the header file. Example usage:
 *
 * @code
 * struct timespec start, stop;
 * Timer(start);
 * // do some work //
 * Timer(stop);
 * double duration = gensvm_elapsed_time(&start, &stop);
 * @endcode
 *
 * This approach to measuring time has been chosen since CLOCK_MONOTONIC_RAW 
 * is guaranteed to be monotonically increasing, and is not affected by leap 
 * seconds or other changes to the system clock. It is therefore thread-safe 
 * and can be used to accurately measure time intervals.
 *
 * @param[in] 	start 	starting time
 * @param[in] 	stop 	end time
 * @returns 		time elapsed in seconds
 */
double gensvm_elapsed_time(struct timespec *start, struct timespec *stop)
{
	double delta_s = stop->tv_sec - start->tv_sec;
	double delta_ns = stop->tv_nsec - start->tv_nsec;
	return delta_s + delta_ns * 1e-9;
}
