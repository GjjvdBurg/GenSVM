/**
 * @file gensvm_print.c
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Various print functions for printing to output streams
 *
 * @details
 * This file contains several utility functions for coordinating input and
 * output of data and model files. It also contains string functions.
 *
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
void err(const char *fmt, ...)
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
