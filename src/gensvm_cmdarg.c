/**
 * @file gensvm_cmdarg.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Functions for dealing with command line arguments
 *
 * @details
 * This file contains several utility functions for coordinating input and
 * output of data and model files.
 *
 */

#include "gensvm_cmdarg.h"

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
