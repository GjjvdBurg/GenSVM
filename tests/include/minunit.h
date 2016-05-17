#undef NDEBUG
#ifndef _minunit_h
#define _minunit_h

#include "dbg.h"
#include <stdlib.h>

#define mu_suite_start() char *message = NULL

#define mu_assert(test, message) if (!(test)) { log_err(message); return message; }

#define mu_run_test(test) debug("\n-----%s", " " #test); \
	message = test(); tests_run++; if (message) return message;

#define RUN_TESTS(name) int main(int argc, char *argv[]) {\
	argc = 1; \
	debug("\n-----\nRUNNING: %s", argv[0]);\
	printf("----\nRUNNING: %s\n", argv[0]);\
	char *result = name();\
	if (result != 0) {\
		printf("\033[91mFAILED:\033[0m %s\n", result);\
	}\
	else {\
		printf("ALL TESTS \033[92mPASSED\033[0m\n");\
	}\
	printf("Tests run: %d\n", tests_run);\
	exit(result != 0);\
}

#define mu_test_missing() printf("\033[33;1mWARNING: Test missing\033[0m\n");

int tests_run;

#endif
