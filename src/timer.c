/**
 * @file timer.c
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Utility functions relating to time
 *
 * @details
 * This file contains a simple function for calculating the time in seconds
 * elapsed between two clock() calls. It also contains a function for
 * generating a string of the current time, used in writing output files.
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

/**
 * @brief Get time string with UTC offset
 *
 * @details
 * Create a string for the current system time. Include an offset of UTC for
 * consistency. The format of the generated string is "DDD MMM D HH:MM:SS 
 * YYYY (UTC +HH:MM)", e.g. "Fri Aug 9, 12:34:56 2013 (UTC +02:00)".
 *
 * @param[in,out] 	buffer 	allocated string buffer, on exit contains
 * 				formatted string
 *
 */
void get_time_string(char *buffer)
{
	int diff, hours, minutes;
	char timestr[MAX_LINE_LENGTH];
	time_t current_time, lt, gt;
	struct tm *lclt;

	// get current time (in epoch)
	current_time = time(NULL);
	if (current_time == ((time_t)-1)) {
		fprintf(stderr, "Failed to compute the current time.\n");
		return;
	}

	// convert time to local time and create a string
	lclt = localtime(&current_time);
	strftime(timestr, MAX_LINE_LENGTH, "%c", lclt);
	if (timestr == NULL) {
		fprintf(stderr, "Failed to convert time to string.\n");
		return;
	}

	// calculate the UTC offset including DST
	lt = mktime(localtime(&current_time));
	gt = mktime(gmtime(&current_time));
	diff = -difftime(gt, lt);
	hours = (diff/3600);
	minutes = (diff%3600)/60;
	if (lclt->tm_isdst == 1)
		hours++;

	sprintf(buffer, "%s (UTC %+03i:%02i)", timestr, hours, minutes);
}
