CC=gcc
CFLAGS=-Wall -Wno-unused-result -Wsign-compare -g -rdynamic -DNDEBUG $(OPTFLAGS)
INCLUDE=-I../include/ -I./include
LIB=-L../lib
LDFLAGS+=-lcblas -llapack -lm -lgensvm

ifneq ($(strip $(shell ldconfig -p | grep libopenblas)),)
override LDFLAGS+=-lopenblas
else ifneq ($(shell ldconfig -p | grep libatlas),)
override LDFLAGS+=-latlas
else
$(error No OpenBLAS or ATLAS found, please install either or alter this Makefile)
endif

TEST_SRC=$(wildcard src/test_*.c)
TESTS=$(patsubst src/%.c,bin/%,$(TEST_SRC))

.PHONY: all

all: $(TESTS)
	bash ./runtests.sh

bin/%: src/%.c
	@echo $<
	@$(CC) $< -o $@ $(CFLAGS) $(INCLUDE) $(LIB) $(LDFLAGS)

valgrind:
	VALGRIND="valgrind --error-limit=no --leak-check=full \
		 --log-file=/tmp/valgrind-%p.log --track-origins=yes \
		 --show-leak-kinds=all -v" $(MAKE)

cover: CFLAGS += --coverage
cover: LDFLAGS += -lgcov
cover: all

clean:
	rm -rf $(TESTS)
	rm -f ./tests.log
	rm -f /tmp/valgrind*.log
	rm -f *.gcno *.gcda *.gcov
