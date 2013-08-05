VERSION=0.1
CC=gcc
CFLAGS=-Wall -O2 -DVERSION=$(VERSION) -g
INCLUDE= -Iinclude/
EXECS=trainMSVMMaj

.PHONY: all clean tar

all: $(EXECS)

override LDFLAGS+=-lblas -llapack -lm

trainMSVMMaj: src/trainMSVMMaj.c src/libMSVMMaj.o src/util.o
	$(CC) -o trainMSVMMaj src/trainMSVMMaj.c src/libMSVMMaj.o src/util.o $(CFLAGS) $(INCLUDE) $(LDFLAGS)

src/libMSVMMaj.o:
	$(CC) -c -o src/libMSVMMaj.o src/libMSVMMaj.c $(CFLAGS) $(INCLUDE)

src/util.o:
	$(CC) -c -o src/util.o src/util.c $(CFLAGS) $(INCLUDE)

clean:
	rm -rf $(EXECS) *.o src/*.o
