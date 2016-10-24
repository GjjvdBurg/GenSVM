/**
 * @file gensvm_timer.c
 * @author G.J.J. van den Burg
 * @date 2014-01-07
 * @brief Utility functions relating to time
 *
 * @details
 * This file contains a simple function for calculating the time in seconds
 * elapsed between two Timer() calls.
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
