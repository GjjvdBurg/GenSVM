VERSION=0.1
CC=gcc
CFLAGS=-Wall -O2 -DVERSION=$(VERSION) -g
INCLUDE= -Iinclude/
EXECS=trainMSVMMaj predMSVMMaj

.PHONY: all clean tar

all: $(EXECS)

override LDFLAGS+=-lblas -llapack -lm

trainMSVMMaj: src/trainMSVMMaj.c src/libMSVMMaj.o src/util.o src/matrix.o
	$(CC) -o trainMSVMMaj src/trainMSVMMaj.c src/libMSVMMaj.o src/util.o src/matrix.o $(CFLAGS) $(INCLUDE) $(LDFLAGS)

predMSVMMaj: src/predMSVMMaj.c src/libMSVMMaj.o src/util.o src/matrix.o
	$(CC) -o predMSVMMaj src/predMSVMMaj.c src/libMSVMMaj.o src/util.o src/matrix.o $(CFLAGS) $(INCLUDE) $(LDFLAGS)

src/libMSVMMaj.o:
	$(CC) -c -o src/libMSVMMaj.o src/libMSVMMaj.c $(CFLAGS) $(INCLUDE)

src/util.o:
	$(CC) -c -o src/util.o src/util.c $(CFLAGS) $(INCLUDE)

src/matrix.o:
	$(CC) -c -o src/matrix.o src/matrix.c $(CFLAGS) $(INCLUDE)

clean:
	rm -rf $(EXECS) *.o src/*.o
