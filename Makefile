VERSION=0.1
CC=gcc
CFLAGS=-Wall -O3 -DVERSION=$(VERSION)
INCLUDE= -Iinclude
LIB= -Llib
EXECS=trainGenSVM trainGenSVMdataset

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

trainGenSVM: src/trainGenSVM.c lib/libgensvm.a
	@$(CC) -o trainGenSVM src/trainGenSVM.c $(CFLAGS) $(INCLUDE) $(LIB)\
		-lgensvm $(LDFLAGS)
	@echo trainGenSVM...

trainGenSVMdataset: src/trainGenSVMdataset.c lib/libgensvm.a
	@$(CC) -o trainGenSVMdataset src/trainGenSVMdataset.c $(CFLAGS) $(INCLUDE) $(LIB) -lgensvm $(LDFLAGS)
	@echo trainGenSVMdataset...

predGenSVM: src/predGenSVM.c lib/libgensvm.a
	@$(CC) -o predGenSVM src/predGenSVM.c $(CFLAGS) $(INCLUDE) $(LIB) -lgensvm $(LDFLAGS)
	@echo predGenSVM...

src/crossval.o:
	@$(CC) -c -o src/crossval.o src/crossval.c $(CFLAGS) $(INCLUDE)
	@echo crossval.o...

src/gensvm_kernel.o:
	@$(CC) -c -o src/gensvm_kernel.o src/gensvm_kernel.c $(CFLAGS) $(INCLUDE)
	@echo gensvm_kernel.o...

src/libGenSVM.o:
	@$(CC) -c -o src/libGenSVM.o src/libGenSVM.c $(CFLAGS) $(INCLUDE)
	@echo libGenSVM.o...

src/gensvm_matrix.o:
	@$(CC) -c -o src/gensvm_matrix.o src/gensvm_matrix.c $(CFLAGS) $(INCLUDE)
	@echo gensvm_matrix.o...

src/gensvm_init.o:
	@$(CC) -c -o src/gensvm_init.o src/gensvm_init.c $(CFLAGS) $(INCLUDE)
	@echo gensvm_init.o...

src/gensvm_io.o:
	@$(CC) -c -o $@ src/gensvm_io.c $(CFLAGS) $(INCLUDE)

src/gensvm_pred.o:
	@$(CC) -c -o src/gensvm_pred.o src/gensvm_pred.c $(CFLAGS) $(INCLUDE)
	@echo gensvm_pred.o...

src/gensvm_sv.o:
	@$(CC) -c -o src/gensvm_sv.o src/gensvm_sv.c $(CFLAGS) $(INCLUDE)
	@echo gensvm_sv.o...

src/gensvm_train.o:
	@$(CC) -c -o src/gensvm_train.o src/gensvm_train.c $(CFLAGS) $(INCLUDE)
	@echo gensvm_train.o...

src/gensvm_train_dataset.o:
	@$(CC) -c -o src/gensvm_train_dataset.o src/gensvm_train_dataset.c $(CFLAGS) $(INCLUDE) 
	@echo gensvm_train_dataset.o...

src/gensvm_lapack.o:
	@$(CC) -c -o src/gensvm_lapack.o src/gensvm_lapack.c $(CFLAGS) $(INCLUDE)
	@echo mylapack.o...

src/strutil.o:
	@$(CC) -c -o src/strutil.o src/strutil.c $(CFLAGS) $(INCLUDE)
	@echo strutil.o...

src/timer.o:
	@$(CC) -c -o src/timer.o src/timer.c $(CFLAGS) $(INCLUDE)
	@echo timer.o...

src/util.o:
	@$(CC) -c -o src/util.o src/util.c $(CFLAGS) $(INCLUDE)
	@echo util.o...

clean:
	rm -rf $(EXECS) *.o src/*.o lib/*.a
