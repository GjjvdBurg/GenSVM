/**
 * @file util.c
 * @author Gertjan van den Burg
 * @date January, 2014
 * @brief Utility functions
 *
 * @details
 * This file contains several utility functions for coordinating input and
 * output of data and model files. It also contains string functions.
 *
 */
#include <stdarg.h>

#include "gensvm_util.h"

FILE *GENSVM_OUTPUT_FILE; 	///< The #GENSVM_OUTPUT_FILE specifies the
				///< output stream to which all output is
				///< written. This is done through the
				///< internal (!)
				///< function gensvm_print_string(). The
				///< advantage of using a global output
				///< stream variable is that the output can
				///< temporarily be suppressed by importing
				///< this variable through @c extern and
				///< (temporarily) setting it to NULL.

/**
 * @brief Check if any command line arguments contain string
 *
 * @details
 * Check if any of a given array of command line arguments contains a given
 * string. If the string is found, the index of the string in argv is
 * returned. If the string is not found, 0 is returned.
 *
 * This function is copied from MSVMpack/libMSVM.c.
 *
 * @param[in] 	argc 	number of command line arguments
 * @param[in] 	argv 	command line arguments
 * @param[in] 	str 	string to find in the arguments
 * @returns 		index of the string in the arguments if found, 0
 * 			otherwise
 */
int gensvm_check_argv(int argc, char **argv, char *str)
{
	int i;
	int arg_str = 0;
	for (i=1; i<argc; i++)
		if (strstr(argv[i], str) != NULL) {
			arg_str = i;
			break;
		}

	return arg_str;
}

/**
 * @brief Check if a command line argument equals a string
 *
 * @details
 * Check if any of the command line arguments is exactly equal to a given
 * string. If so, return the index of the corresponding command line argument.
 * If not, return 0.
 *
 * This function is copied from MSVMpack/libMSVM.c
 *
 * @param[in] 	argc 	number of command line arguments
 * @param[in] 	argv 	command line arguments
 * @param[in] 	str 	string to find in the arguments
 * @returns 		index of the command line argument that corresponds to
 * 			the string, 0 if none matches.
 */
int gensvm_check_argv_eq(int argc, char **argv, char *str)
{
	int i;
	int arg_str = 0;
	for (i=1; i<argc; i++)
		if (strcmp(argv[i], str) == 0) {
			arg_str = i;
			break;
		}

	return arg_str;
}


/**
 * @brief Print a given string to the specified output stream
 *
 * @details
 * This function is used to print a given string to the output stream
 * specified by #GENSVM_OUTPUT_FILE. The stream is flushed after the string
 * is written to the stream. If #GENSVM_OUTPUT_FILE is NULL, nothing is
 * written. Note that this function is only used by note(), it should never be
 * used directly.
 *
 * @param[in] 	s 	string to write to the stream
 *
 */
static void gensvm_print_string(const char *s)
{
	if (GENSVM_OUTPUT_FILE != NULL) {
		fputs(s, GENSVM_OUTPUT_FILE);
		fflush(GENSVM_OUTPUT_FILE);
	}
}

/**
 * @brief Parse a formatted string and write to the output stream
 *
 * @details
 * This function is a replacement of fprintf(), such that the output stream
 * does not have to be specified at each function call. The functionality is
 * exactly the same however. Writing the formatted string to the output stream
 * is handled by gensvm_print_string().
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
	(*gensvm_print_string)(buf);
}
