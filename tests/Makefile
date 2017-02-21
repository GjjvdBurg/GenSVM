CC=gcc
CFLAGS=-Wall -g -rdynamic -DNDEBUG $(OPTFLAGS)
INCLUDE=-I../include/ -I./include
LIB=-L../lib

ifdef NOATLAS
override LDFLAGS+=-lm -lcblas -llapack -lgensvm
else
override LDFLAGS+=-lm -lcblas -llapack -latlas -lgensvm
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
