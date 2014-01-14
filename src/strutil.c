/**
 * @file strutil.c
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Utility functions for dealing with strings
 *
 * @details
 * This file contains functions for reading files, reading strings from a
 * format and checking start and ends of strings.
 */

#include "strutil.h"

/**
 * @brief Check if a string starts with a prefix
 *
 * @param[in] 	str 	string
 * @param[in] 	pre 	prefix
 * @returns 		boolean, true if string starts with prefix, false
 * 			otherwise
 */
bool str_startswith(const char *str, const char *pre)
{
	size_t lenpre = strlen(pre),
	       lenstr = strlen(str);
	return lenstr < lenpre ? false : strncmp(pre, str, lenpre) == 0;
}

/**
 * @brief Check if a string ends with a suffix
 *
 * @param[in] 	str 	string
 * @param[in] 	suf 	suffix
 * @returns 		boolean, true if string ends with suffix, false
 * 			otherwise
 */
bool str_endswith(const char *str, const char *suf)
{
	size_t lensuf = strlen(suf),
	       lenstr = strlen(str);
	return lenstr < lensuf ? false : strncmp(str + lenstr - lensuf, suf, 
			lensuf) == 0;
}

/**
 * @brief Move to next line in file
 *
 * @param[in] 	fid 		File opened for reading
 * @param[in] 	filename 	name of the file pointed to by fid
 */
void next_line(FILE *fid, char *filename)
{
	char buffer[MAX_LINE_LENGTH];
	get_line(fid, filename, buffer);
}

/**
 * @brief Read line to buffer
 *
 * @param[in] 		fid 		File opened for reading
 * @param[in] 		filename 	name of the file
 * @param[in,out] 	buffer 		allocated buffer to read to
 */
void get_line(FILE *fid, char *filename, char *buffer)
{
	if (fgets(buffer, MAX_LINE_LENGTH, fid) == NULL) {
		fprintf(stderr, "Error reading file %s\n", filename);
		exit(1);
	}
}

/**
 * @brief Read a double from file following a format
 *
 * @param[in] 	fid 		File opened for reading
 * @param[in] 	filename 	Name of the file
 * @param[in] 	fmt 		Format containing a float format
 * @returns 			value read (if any)
 */
double get_fmt_double(FILE *fid, char *filename, const char *fmt)
{
	char buffer[MAX_LINE_LENGTH];
	double value;

	get_line(fid, filename, buffer);
	sscanf(buffer, fmt, &value);
	return value;
}

/**
 * @brief Read a long integer from file following a format
 *
 * @param[in] 	fid 		File opened for reading
 * @param[in] 	filename 	Name of the file
 * @param[in] 	fmt 		Format containing a long integer format
 * @returns 			value read (if any)
 */
long get_fmt_long(FILE *fid, char *filename, const char *fmt)
{
	char buffer[MAX_LINE_LENGTH];
	long value;

	get_line(fid, filename, buffer);
	sscanf(buffer, fmt, &value);
	return value;
}

/**
 * @brief Read all doubles in a given buffer
 *
 * @details
 * This function is used to read a line of doubles from a buffer. All the
 * doubles found are stored in a pre-allocated array. 
 *
 * @param[in] 	buffer 		a string buffer
 * @param[in] 	offset 		an offset of the string to start looking for
 * 				doubles
 * @param[in] 	all_doubles 	pre-allocated array of doubles (should be large
 * 				enough)
 * @returns 			number of doubles read
 */
long all_doubles_str(char *buffer, long offset, double *all_doubles)
{
	double value;
	long i = 0;
	char *start, *end;

	start = buffer + offset;
	while (true) {
		value = strtod(start, &end);
		if (start != end) {
			all_doubles[i] = value;
			i++;
		} else 
			break;
		start = end;
		end = NULL;
	}

	return i;
}

/**
 * @brief Read all longs in a given buffer
 *
 * @details
 * This function is used to read a line of longs from a buffer. All the
 * longs found are stored in a pre-allocated array. 
 *
 * @param[in] 	buffer 		a string buffer
 * @param[in] 	offset 		an offset of the string to start looking for
 * 				longs
 * @param[in] 	all_longs 	pre-allocated array of longs (should be large
 * 				enough)
 * @returns 			number of longs read
 */
long all_longs_str(char *buffer, long offset, long *all_longs)
{
	long value;
	long i = 0;
	char *start, *end;

	start = buffer + offset;
	while (true) {
		value = strtol(start, &end, 10);
		if (start != end) {
			all_longs[i] = value;
			i++;
		} else
			break;
		start = end;
		end = NULL;
	}

	return i;
}
