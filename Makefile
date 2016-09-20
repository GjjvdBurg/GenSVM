VERSION=0.1
CC=gcc
CFLAGS=-Wall -DVERSION=$(VERSION) -g -O3
INCLUDE= -Iinclude
LIB= -Llib
DOXY=doxygen
DOCDIR=doc
DOXYFILE=$(DOCDIR)/Doxyfile
LCOV=lcov
GENHTML=genhtml

EXECS=gensvm gensvm_grid

.PHONY: all clean doc test cover

all: lib/libgensvm.a $(EXECS)

override LDFLAGS+=-lcblas -llapack -lm

debug: CFLAGS += -DDEBUG
debug: all

doc:
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

lib/libgensvm.a: \
	src/gensvm_base.o \
	src/gensvm_cmdarg.o \
	src/gensvm_copy.o \
	src/gensvm_cv_util.o \
	src/gensvm_debug.o \
	src/gensvm_grid.o \
	src/gensvm_gridsearch.o \
	src/gensvm_init.o \
	src/gensvm_io.o \
	src/gensvm_kernel.o \
	src/gensvm_memory.o \
	src/gensvm_optimize.o \
	src/gensvm_pred.o \
	src/gensvm_print.o \
	src/gensvm_queue.o \
	src/gensvm_simplex.o \
	src/gensvm_strutil.o \
	src/gensvm_sv.o \
	src/gensvm_task.o \
	src/gensvm_timer.o \
	src/gensvm_train.o
	@ar rcs lib/libgensvm.a \
		src/gensvm_base.o \
		src/gensvm_cmdarg.o \
		src/gensvm_copy.o \
		src/gensvm_cv_util.o \
		src/gensvm_debug.o \
		src/gensvm_grid.o \
		src/gensvm_gridsearch.o \
		src/gensvm_init.o \
		src/gensvm_io.o \
		src/gensvm_kernel.o \
		src/gensvm_memory.o \
		src/gensvm_optimize.o \
		src/gensvm_pred.o \
		src/gensvm_print.o \
		src/gensvm_queue.o \
		src/gensvm_simplex.o \
		src/gensvm_strutil.o \
		src/gensvm_sv.o \
		src/gensvm_task.o \
		src/gensvm_timer.o \
		src/gensvm_train.o
	@echo libgensvm.a...

gensvm: src/GenSVMtraintest.c lib/libgensvm.a
	@$(CC) -o $@ $< $(CFLAGS) $(INCLUDE) $(LIB) -lgensvm $(LDFLAGS)
	@echo gensvm ...

gensvm_grid: src/GenSVMgrid.c lib/libgensvm.a
	@$(CC) -o $@ $< $(CFLAGS) $(INCLUDE) $(LIB) -lgensvm $(LDFLAGS)
	@echo gensvm_grid ...

src/%.o: src/%.c
	@$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -c $< -o $@
	@echo $< ...
