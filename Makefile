VERSION=0.1
CC=gcc
CFLAGS=-Wall -O3 -DVERSION=$(VERSION)
INCLUDE= -Iinclude
LIB= -Llib

EXECS=GenSVM_train GenSVM_grid gensvm

.PHONY: all clean tar

all: lib/libgensvm.a $(EXECS)

override LDFLAGS+=-lcblas -llapack -lm

lib/libgensvm.a: \
	src/crossval.o \
	src/libGenSVM.o \
	src/gensvm_init.o \
	src/gensvm_io.o \
	src/gensvm_kernel.o \
	src/gensvm_lapack.o \
	src/gensvm_matrix.o \
	src/gensvm_pred.o \
	src/gensvm_sv.o \
	src/gensvm_train.o \
	src/gensvm_train_dataset.o \
	src/strutil.o \
	src/timer.o \
	src/util.o 
	@ar rcs lib/libgensvm.a \
		src/crossval.o \
		src/libGenSVM.o \
		src/gensvm_init.o \
		src/gensvm_io.o \
		src/gensvm_matrix.o \
		src/gensvm_kernel.o \
		src/gensvm_lapack.o \
		src/gensvm_pred.o \
		src/gensvm_sv.o \
		src/gensvm_train.o \
		src/gensvm_train_dataset.o \
		src/strutil.o \
		src/timer.o \
		src/util.o 
	@echo libgensvm.a...

gensvm: src/GenSVMtraintest.c lib/libgensvm.a
	@$(CC) -o $@ $< $(CFLAGS) $(INCLUDE) $(LIB) -lgensvm $(LDFLAGS)
	@echo gensvm...

GenSVM_train: src/GenSVMtrain.c lib/libgensvm.a
	@$(CC) -o GenSVM_train src/GenSVMtrain.c $(CFLAGS) $(INCLUDE) $(LIB)\
		-lgensvm $(LDFLAGS)
	@echo GenSVM_train...

GenSVM_grid: src/GenSVMgrid.c lib/libgensvm.a
	@$(CC) -o GenSVM_grid src/GenSVMgrid.c $(CFLAGS) $(INCLUDE) $(LIB) \
		-lgensvm $(LDFLAGS)
	@echo GenSVM_grid...

GenSVM_pred: src/GenSVMpred.c lib/libgensvm.a
	@$(CC) -o GenSVM_pred src/GenSVMpred.c $(CFLAGS) $(INCLUDE) $(LIB) \
		-lgensvm $(LDFLAGS)
	@echo GenSVM_pred...

clean:
	rm -rf $(EXECS) *.o src/*.o lib/*.a

src/%.o: src/%.c
	@$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -c $< -o $@
	@echo $<...
