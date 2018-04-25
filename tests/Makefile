CC=gcc
CFLAGS=-Wall -Wno-unused-result -Wsign-compare -g -rdynamic -DNDEBUG $(OPTFLAGS) -fsanitize=undefined
INCLUDE=-I../include/ -I./include
LIB=-L../lib

#############################################
#    COLUMN-MAJOR ORDER MUST BE SET HERE    #
#  (should be the same as for the library)  #
#############################################
COLUMN_MAJOR=false
ifeq ($(COLUMN_MAJOR), true)
	CFLAGS += -DCOLUMN_MAJOR_ORDER
endif

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
