VERSION=0.1
CC=gcc
CFLAGS=-Wall -O3 -DVERSION=$(VERSION) -g
INCLUDE= -Iinclude
LIB= -Llib
EXECS=trainMSVMMaj trainMSVMMajdataset

.PHONY: all clean tar

all: lib/libmsvmmaj.a $(EXECS)

override LDFLAGS+=-lblas -llapack -lm -lcblas

lib/libmsvmmaj.a: \
	src/crossval.o \
	src/libMSVMMaj.o \
	src/msvmmaj_init.o \
	src/msvmmaj_io.o \
	src/msvmmaj_kernel.o \
	src/msvmmaj_lapack.o \
	src/msvmmaj_matrix.o \
	src/msvmmaj_pred.o \
	src/msvmmaj_train.o \
	src/msvmmaj_train_dataset.o \
	src/strutil.o \
	src/timer.o \
	src/util.o 
	@ar rcs lib/libmsvmmaj.a \
		src/crossval.o \
		src/libMSVMMaj.o \
		src/msvmmaj_init.o \
		src/msvmmaj_io.o \
		src/msvmmaj_matrix.o \
		src/msvmmaj_kernel.o \
		src/msvmmaj_lapack.o \
		src/msvmmaj_pred.o \
		src/msvmmaj_train.o \
		src/msvmmaj_train_dataset.o \
		src/strutil.o \
		src/timer.o \
		src/util.o 
	@echo libmsvmmaj.a...

trainMSVMMaj: src/trainMSVMMaj.c lib/libmsvmmaj.a
	@$(CC) -o trainMSVMMaj src/trainMSVMMaj.c $(CFLAGS) $(INCLUDE) $(LIB) -lmsvmmaj $(LDFLAGS)
	@echo trainMSVMMaj...

trainMSVMMajdataset: src/trainMSVMMajdataset.c lib/libmsvmmaj.a
	@$(CC) -o trainMSVMMajdataset src/trainMSVMMajdataset.c $(CFLAGS) $(INCLUDE) $(LIB) -lmsvmmaj $(LDFLAGS)
	@echo trainMSVMMajdataset...

predMSVMMaj: src/predMSVMMaj.c lib/libmsvmmaj.a
	@$(CC) -o predMSVMMaj src/predMSVMMaj.c $(CFLAGS) $(INCLUDE) $(LIB) -lmsvmmaj $(LDFLAGS)
	@echo predMSVMMaj...

src/crossval.o:
	@$(CC) -c -o src/crossval.o src/crossval.c $(CFLAGS) $(INCLUDE)
	@echo crossval.o...

src/msvmmaj_kernel.o:
	@$(CC) -c -o src/msvmmaj_kernel.o src/msvmmaj_kernel.c $(CFLAGS) $(INCLUDE)
	@echo msvmmaj_kernel.o...

src/libMSVMMaj.o:
	@$(CC) -c -o src/libMSVMMaj.o src/libMSVMMaj.c $(CFLAGS) $(INCLUDE)
	@echo libMSVMMaj.o...

src/msvmmaj_matrix.o:
	@$(CC) -c -o src/msvmmaj_matrix.o src/msvmmaj_matrix.c $(CFLAGS) $(INCLUDE)
	@echo msvmmaj_matrix.o...

src/msvmmaj_init.o:
	@$(CC) -c -o src/msvmmaj_init.o src/msvmmaj_init.c $(CFLAGS) $(INCLUDE)
	@echo msvmmaj_init.o...

src/msvmmaj_io.o:
	@$(CC) -c -o $@ src/msvmmaj_io.c $(CFLAGS) $(INCLUDE)

src/msvmmaj_pred.o:
	@$(CC) -c -o src/msvmmaj_pred.o src/msvmmaj_pred.c $(CFLAGS) $(INCLUDE)
	@echo msvmmaj_pred.o...

src/msvmmaj_train.o:
	@$(CC) -c -o src/msvmmaj_train.o src/msvmmaj_train.c $(CFLAGS) $(INCLUDE)
	@echo msvmmaj_train.o...

src/msvmmaj_train_dataset.o:
	@$(CC) -c -o src/msvmmaj_train_dataset.o src/msvmmaj_train_dataset.c $(CFLAGS) $(INCLUDE) 
	@echo msvmmaj_train_dataset.o...

src/msvmmaj_lapack.o:
	@$(CC) -c -o src/msvmmaj_lapack.o src/msvmmaj_lapack.c $(CFLAGS) $(INCLUDE)
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
