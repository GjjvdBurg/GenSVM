VERSION=0.1
CC=gcc
CFLAGS=-Wall -O3 -DVERSION=$(VERSION) -g
INCLUDE= -Iinclude
LIB= -Llib
DOXY=doxygen
DOCDIR=doc
DOXYFILE=$(DOCDIR)/Doxyfile

EXECS=gensvm gensvm_grid

.PHONY: all clean doc test

all: lib/libgensvm.a $(EXECS)

override LDFLAGS+=-lcblas -llapack -lm

debug: CFLAGS += -DDEBUG
debug: all

doc:
	$(DOXY) $(DOXYFILE)

clean:
	rm -rf $(EXECS) *.o src/*.o lib/*.a
	$(MAKE) -C tests clean

test: lib/libgensvm.a
	$(MAKE) -C tests all

lib/libgensvm.a: \
	src/gensvm_base.o \
	src/gensvm_cmdarg.o \
	src/gensvm_copy.o \
	src/gensvm_cv_util.o \
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
	src/gensvm_timer.o
	@ar rcs lib/libgensvm.a \
		src/gensvm_base.o \
		src/gensvm_cmdarg.o \
		src/gensvm_copy.o \
		src/gensvm_cv_util.o \
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
		src/gensvm_timer.o
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
