# Makefile for GenSVM
#
# Author: G.J.J. van den Burg
# Copyright (c) 2016 G.J.J. van den Burg
# License: GPLv2
#

SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules

VERSION=0.2.2
CC=gcc
CFLAGS=-Wall -Wno-unused-result -Wsign-compare -Wstrict-prototypes \
       -DVERSION=$(VERSION) -g -O3
INCLUDE= -Iinclude
LIB= -Llib
DOXY=doxygen
DOCDIR=doc
DOXYFILE=$(DOCDIR)/Doxyfile
LCOV=lcov
GENHTML=genhtml
LDFLAGS+=-lcblas -llapack -lm

EXECS=gensvm gensvm_grid

# Should be a cleaner way to do this if we rename the exec sources
EXECS_C=src/GenSVMtraintest.c src/GenSVMgrid.c
SRC=$(filter-out $(EXECS_C),$(wildcard src/*.c))
OBJ=$(patsubst %.c,%.o,$(SRC))

.PHONY: all clean doc test cover

all: lib/libgensvm.a $(EXECS)

ifneq ($(strip $(shell ldconfig -p | grep libopenblas)),)
override LDFLAGS+=-lopenblas
else ifneq ($(shell ldconfig -p | grep libatlas),)
override LDFLAGS+=-latlas
else
$(error No OpenBLAS or ATLAS found, please install either or alter this Makefile. If you believe this is a bug please open an issue at https://www.github.com/GjjvdBurg/GenSVM)
endif

debug: CFLAGS += -DDEBUG
debug: all

doc: cover
	$(DOXY) $(DOXYFILE)

clean:
	rm -rf $(EXECS) *.o src/*.o lib/*.a *.{gcno,gcov} src/*.{gcno,gcda}
	$(MAKE) -C tests clean

test: lib/libgensvm.a
	$(MAKE) -C tests all

cover: CFLAGS += --coverage
cover: LDFLAGS += --coverage -lgcov
cover: CFLAGS := $(filter-out -O3,$(CFLAGS))
cover: lib/libgensvm.a
	$(LCOV) -c -i -d ./src/ -o ./cover/coverage.base
	$(MAKE) -C tests cover
	mkdir -p cover
	$(LCOV) -c -d ./src/ -o ./cover/coverage.run
	$(LCOV) -a ./cover/coverage.base -a ./cover/coverage.run \
		-o ./cover/coverage.all
	$(GENHTML) -o ./cover ./cover/coverage.all
	rm -f src/*.{gcda,gcno} tests/*.{gcda,gcno}
	cp -r cover doc/html/

lib/libgensvm.a: $(OBJ)
	@ar rcs lib/libgensvm.a $(OBJ)
	@echo libgensvm.a...

gensvm: src/GenSVMtraintest.c lib/libgensvm.a
	$(CC) -o $@ $< $(CFLAGS) $(INCLUDE) $(LIB) -lgensvm $(LDFLAGS)

gensvm_grid: src/GenSVMgrid.c lib/libgensvm.a
	$(CC) -o $@ $< $(CFLAGS) $(INCLUDE) $(LIB) -lgensvm $(LDFLAGS)

src/%.o: src/%.c
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -c $< -o $@
