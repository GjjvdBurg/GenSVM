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
 * @brief Cross-platform time recording
 *
 * @details
 * This function computes the current time on both Windows and Unix-like 
 * platforms. On Windows the QueryPerformanceCounter() function is used and on 
 * unix the gettimeofday() function.
 *
 * @param[out] 	t 	a GenTime type to record the time in
 *
 */
void gensvm_set_time(GenTime *t)
{
	#ifdef ON_WINDOWS
	if (!QueryPerformanceCounter(t))
		FatalError("QueryPerformanceCounter failed.");
	#else
	gettimeofday(t, 0);
	#endif
}

/**
 * @brief Calculate the time between two time recordings
 *
 * @details
 * This function should be used with time recordings done by 
 * gensvm_set_time(). To do this easily the Timer() macro has been defined in 
 * the header file.  Example usage:
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
double gensvm_elapsed_time(GenTime *start, GenTime *stop)
{
	#ifdef ON_WINDOWS
	static LARGE_INTEGER freq;
	if (!freq.QuadPart) {
		if (!QueryPerformanceFrequency(&freq))
			FatalError("QueryPerformanceFrequency failed.");
	}
	return (double)(stop->QuadPart - start->QuadPart) / freq.QuadPart;
	#else
	double diff_sec = stop->tv_sec - start->tv_sec;
	double diff_usec = stop->tv_usec - start->tv_usec;
	return diff_sec + diff_usec / 1000000.;
	#endif
}


// https://stackoverflow.com/a/23159148
/**
 * @brief Cross-platform sleep function
 *
 * @param[in] 	seconds 	time to sleep in seconds. This can be a double
 * 				but will be rounded to the nearest
 * 				millisecond.
 */
void gensvm_sleep(double seconds)
{
	int milliseconds = 1000.0 * seconds;
	#ifdef ON_WINDOWS
	Sleep(milliseconds);
	#else
	usleep(1000 * milliseconds);
	#endif
}
