/**
 * @file gensvm_strutil.c
 * @author G.J.J. van den Burg
 * @date 2014-01-07
 * @brief Utility functions for dealing with strings
 *
 * @details
 * This file contains functions for reading files, reading strings from a
 * format and checking start and ends of strings.
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

#include "gensvm_strutil.h"
#include "gensvm_print.h"

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
 * @brief Check if a string contains a char
 *
 * @details
 * Simple utility function to check if a char occurs in a string.
 *
 * @param[in] 	str 	input string
 * @param[in] 	c 	character
 *
 * @return 		number of times c occurs in str
 */
bool str_contains_char(const char *str, const char c)
{
	size_t i, len = strlen(str);
	for (i=0; i<len; i++)
		if (str[i] == c)
			return true;
	return false;
}

/**
 * @brief Count the number of times a string contains any character of another
 *
 * @details
 * This function is used to count the number of expected parts in the function 
 * str_split(). It counts the number of times a character from a string of 
 * characters is present in an input string.
 *
 * @param[in] 	str 	input string
 * @param[in] 	chars 	characters to count
 *
 * @return  		number of times any character from chars occurs in str
 *
 */
int count_str_occurrences(const char *str, const char *chars)
{
	size_t i, j, len_str = strlen(str),
	    len_chars = strlen(chars);
	int count = 0;
	for (i=0; i<len_str; i++) {
		for (j=0; j<len_chars; j++) {
			count += (str[i] == chars[j]);
		}
	}
	return count;
}

/**
 * @brief Split a string on delimiters and return an array of parts
 *
 * @details
 * This function takes as input a string and a string of delimiters. As 
 * output, it gives an array of the parts of the first string, splitted on the 
 * characters in the second string. The input string is not changed, and the 
 * output contains all copies of the input string parts.
 *
 * @note
 * The code is based on: http://stackoverflow.com/a/9210560
 *
 * @param[in] 	original 	string you wish to split
 * @param[in] 	delims 		string with delimiters to split on
 * @param[out] 	len_ret 	length of the output array
 *
 * @return 			array of string parts
 */
char **str_split(char *original, const char *delims, int *len_ret)
{
	char *copy = NULL,
	     *token = NULL,
	     **result = NULL;
	size_t i, count;
	size_t len = strlen(original);
	size_t n_delim = strlen(delims);

	// add the null terminator to the delimiters
	char all_delim[1 + n_delim];
	for (i=0; i<n_delim; i++)
		all_delim[i] = delims[i];
	all_delim[n_delim] = '\0';

	// number of occurances of the delimiters
	count = count_str_occurrences(original, delims);

	// extra count in case there is a delimiter at the end
	count += (str_contains_char(delims, original[len - 1]));

	// extra count for the null terminator
	count++;

	// allocate the result array
	result = Calloc(char *, count);

	// tokenize a copy of the original string and keep the splits
	i = 0;
	copy = Calloc(char, len + 1);
	strcpy(copy, original);
	token = strtok(copy, all_delim);
	while (token) {
		result[i] = Calloc(char, strlen(token) + 1);
		strcpy(result[i], token);
		i++;

		token = strtok(NULL, all_delim);
	}
	free(copy);

	*len_ret = i;

	return result;
}

/**
 * @brief Move to next line in file
 *
 * @param[in] 	fid 		File opened for reading
 * @param[in] 	filename 	name of the file pointed to by fid
 */
void next_line(FILE *fid, char *filename)
{
	char buffer[GENSVM_MAX_LINE_LENGTH];
	get_line(fid, filename, buffer);
}

/**
 * @brief Read line to buffer
 *
 * @param[in] 		fid 		File opened for reading
 * @param[in] 		filename 	name of the file
 * @param[in,out] 	buffer 		allocated buffer to read to
 */
char *get_line(FILE *fid, char *filename, char *buffer)
{
	char *retval = fgets(buffer, GENSVM_MAX_LINE_LENGTH, fid);
	if (retval == NULL) {
		gensvm_error("[GenSVM Error]: Error reading from file %s\n", filename);
	}
	return retval;
}

/**
 * @brief Read a double from file following a format
 *
 * @details
 * This function reads a double value from a file. If no value can be found, a
 * warning is printed to stderr, and NAN is returned.
 *
 * @param[in] 	fid 		File opened for reading
 * @param[in] 	filename 	Name of the file
 * @param[in] 	fmt 		Format containing a float format
 * @returns 			value read (if any)
 */
double get_fmt_double(FILE *fid, char *filename, const char *fmt)
{
	char buffer[GENSVM_MAX_LINE_LENGTH];
	double value = NAN;
	int retval;

	get_line(fid, filename, buffer);
	retval = sscanf(buffer, fmt, &value);
	if (retval == 0)
		gensvm_error("[GenSVM Error]: No double read from file.\n");
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
	char buffer[GENSVM_MAX_LINE_LENGTH];
	long value = 0;
	int retval;

	get_line(fid, filename, buffer);
	retval = sscanf(buffer, fmt, &value);
	if (retval == 0)
		gensvm_error("[GenSVM Error]: No long read from file.\n");
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
	char *start = NULL,
	     *end = NULL;

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
	char *start = NULL,
	     *end = NULL;

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
