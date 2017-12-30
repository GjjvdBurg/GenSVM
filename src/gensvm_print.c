/**
 * @file gensvm_print.c
 * @author G.J.J. van den Burg
 * @date 2014-01-07
 * @brief Various print functions for printing to output streams
 *
 * @details
 * This file contains several utility functions for coordinating input and
 * output of data and model files. It also contains string functions.
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

#include "gensvm_print.h"

FILE *GENSVM_OUTPUT_FILE = NULL; 	///< The #GENSVM_OUTPUT_FILE specifies the
				///< output stream to which all output is
				///< written. This is done through the
				///< function note(). The
				///< advantage of using a global output
				///< stream variable is that the output can
				///< temporarily be suppressed by importing
				///< this variable through @c extern and
				///< (temporarily) setting it to NULL.

FILE *GENSVM_ERROR_FILE = NULL; 	///< The #GENSVM_ERROR_FILE specifies the
				///< output stream to use when writing an
				///< error.  Typically this is stderr, but
				///< when unit testing we can temporarily
				///< redirect this to check if the correct
				///< output is written.

/**
 * @brief Parse a formatted string and write to the output stream
 *
 * @details
 * This function is a replacement of fprintf(), such that the output stream
 * does not have to be specified at each function call. The functionality is
 * exactly the same however.
 *
 * @param[in] 	fmt 	String format
 * @param[in] 	... 	variable argument list for the string format
 *
 */
void note(const char *fmt,...)
{
	char buf[BUFSIZ];
	va_list ap;
	va_start(ap,fmt);
	vsprintf(buf,fmt,ap);
	va_end(ap);
	if (GENSVM_OUTPUT_FILE != NULL) {
		fputs(buf, GENSVM_OUTPUT_FILE);
		fflush(GENSVM_OUTPUT_FILE);
	}
}

/**
 * @brief Parse a formatted string and write it to standard error
 *
 * @details
 * Shorthand for fprintf(GENSVM_ERROR_FILE, ...)
 *
 * @param[in] 	fmt 	string format
 * @param[in] 	... 	variable argument list for the string format
 */
void gensvm_error(const char *fmt, ...)
{
	char buf[BUFSIZ];
	va_list ap;
	va_start(ap, fmt);
	vsprintf(buf, fmt, ap);
	va_end(ap);
	if (GENSVM_ERROR_FILE != NULL) {
		fputs(buf, GENSVM_ERROR_FILE);
		fflush(GENSVM_ERROR_FILE);
	}
}
