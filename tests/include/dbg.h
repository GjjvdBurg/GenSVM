/**
 * @file dbg.h
 * @brief Debug macros for the minunit framework
 * @author Zed Shaw
 *
 * @details
 * These debug macros come from Zed Shaw's book Learn C The Hard Way, and are
 * used for the testing framework of GenSVM.
 *
 * @sa minunit.h
 */

#ifndef __dbg_h__
#define __dbg_h__

#include <stdio.h>
#include <errno.h>
#include <string.h>

#ifdef NDEBUG
  /**
   * Define debug macro as doing nothing when we don't want debug output
   */
  #define debug(M, ...)
#else
  /**
   * Print debug info to stderr
   */
  #define debug(M, ...) fprintf(stderr, "DEBUG %s:%d: " M "\n", __FILE__, \
  		__LINE__, ##__VA_ARGS__)
#endif

/**
 * Return a clean string of the current errno
 */
#define clean_errno() (errno == 0 ? "None" : strerror(errno))

/**
 * Log an error to stderr
 */
#define log_err(M, ...) fprintf(stderr, "[ERROR] (%s:%d: errno: %s) " M "\n", \
		__FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)

/**
 * Log a warning to stderr
 */
#define log_warn(M, ...) fprintf(stderr, "[WARN] (%s:%d: errno: %s) " M "\n", \
		__FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)

/**
 * Log info to stderr
 */
#define log_info(M, ...) fprintf(stderr, "[INFO] (%s:%d) " M "\n", \
		__FILE__, __LINE__, ##__VA_ARGS__)

/**
 * Check a condition an log an error with the given message if it fails
 */
#define check(A, M, ...) if(!(A)) { log_err(M, ##__VA_ARGS__); errno=0; \
	goto error; }

/**
 * Log an error with the given message and reset errno
 */
#define sentinel(M, ...) { log_err(M, ##__VA_ARGS__); errno=0; goto error; }

/**
 * Check a memory allocation
 */
#define check_mem(A) check((A), "Out of memory.");

/**
 * Check a condition and log to debug if it fails, and reset errno
 */
#define check_debug(A, M, ...) if (!(A)) { debug(M, ##__VA_ARGS__); errno=0; \
  	goto error; }

#endif
