/**
 * @file gensvm_memory.c
 * @author G.J.J. van den Burg
 * @date 2016-05-01
 * @brief Utility functions for memory allocation
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

#include "gensvm_globals.h" // imports gensvm_memory.h

/**
 * @brief Wrapper for calloc() which warns when allocation fails
 *
 * @details
 * This is a wrapper function around calloc from <stdlib.h>. It tries to
 * allocate the requested memory and checks if the memory was correctly
 * allocated. If not, an error is printed using gensvm_error(), which 
 * describes the file and linenumber and size failed to allocate. After this, 
 * the program exits. See also the defines in gensvm_memory.h.
 *
 * @note
 * This function should not be used directly. Calloc() should be used.
 *
 * @param[in] 	file 		filename used for error printing
 * @param[in] 	line 		line number used for error printing
 * @param[in] 	size 		the size to allocate
 * @param[in] 	typesize 	the size of the type to allocate
 * @return 			the pointer to the memory allocated
 */
void *mycalloc(const char *file, int line, unsigned long size,
	       	size_t typesize)
{
	void *ptr = calloc(size, typesize);

	if (!ptr) {
		// LCOV_EXCL_START
		fprintf(stderr, "[GenSVM Error]: Couldn't allocate memory: "
				"%lu bytes (%s:%d)\n", size, file, line);
		exit(EXIT_FAILURE);
		// LCOV_EXCL_STOP
	}
	return ptr;
}

/**
 * @brief Wrapper for malloc() which warns when allocation fails
 *
 * @details
 * This is a wrapper function around malloc from <stdlib.h>. It tries to
 * allocate the requested memory and checks if the memory was correctly
 * allocated. If not, an error is printed using gensvm_error(), which 
 * describes the file and linenumber and size failed to allocate. After this, 
 * the program exits. See also the defines in gensvm_memory.h.
 *
 * @note
 * This function should not be used directly. Malloc() should be used.
 *
 * @param[in] 	file 		filename used for error printing
 * @param[in] 	line 		line number used for error printing
 * @param[in] 	size 		the size to allocate
 * @return 			the pointer to the memory allocated
 */
void *mymalloc(const char *file, int line, unsigned long size)
{
	void *ptr = malloc(size);
	if (!ptr) {
		// LCOV_EXCL_START
		fprintf(stderr, "[GenSVM Error]: Couldn't allocate memory: "
				"%lu bytes (%s:%d)\n", size, file, line);
		exit(EXIT_FAILURE);
		// LCOV_EXCL_STOP
	}
	return ptr;
}

/**
 * @brief Wrapper for realloc() which warns when allocation fails
 *
 * @details
 * This is a wrapper function around realloc from <stdlib.h>. It tries to
 * reallocate the requested memory and checks if the memory was correctly
 * reallocated. If not, an error is printed using gensvm_error(), which 
 * describes the file and linenumber and size failed to reallocate. After 
 * this, the program exits. See also the defines in gensvm_memory.h.
 *
 * @note
 * This function should not be used directly. Realloc() should be used.
 *
 * @param[in] 	file 		filename used for error printing
 * @param[in] 	line 		line number used for error printing
 * @param[in] 	size 		the size to allocate
 * @param[in] 	var 		the pointer to the memory to reallocate
 * @return 			the pointer to the memory reallocated
 */
void *myrealloc(const char *file, int line, unsigned long size, void *var)
{
	void *ptr = realloc(var, size);
	if (!ptr) {
		// LCOV_EXCL_START
		fprintf(stderr, "[GenSVM Error]: Couldn't reallocate memory: "
				"%lu bytes (%s:%d)\n", size, file, line);
		exit(EXIT_FAILURE);
		// LCOV_EXCL_STOP
	}
	return ptr;
}
