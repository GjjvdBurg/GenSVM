/**
 * @file gensvm_memory.c
 * @author Gertjan van den Burg
 * @date May, 2016
 * @brief Utility functions for memory allocation
 *
 */

#include "gensvm_globals.h" // imports gensvm_memory.h

/**
 * @brief Wrapper for calloc() which warns when allocation fails
 *
 * @details
 * This is a wrapper function around calloc from <stdlib.h>. It tries to
 * allocate the requested memory and checks if the memory was correctly
 * allocated. If not, an error is printed using err(), which describes the
 * file and linenumber and size failed to allocate. After this, the program
 * exits. See also the defines in gensvm_memory.h.
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
 * allocated. If not, an error is printed using err(), which describes the
 * file and linenumber and size failed to allocate. After this, the program
 * exits. See also the defines in gensvm_memory.h.
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
 * reallocated. If not, an error is printed using err(), which describes the
 * file and linenumber and size failed to reallocate. After this, the program
 * exits. See also the defines in gensvm_memory.h.
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
