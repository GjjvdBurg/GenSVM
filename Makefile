VERSION=0.2.1
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
EXECS=gensvm gensvm_grid

#############################################
#     COLUMN-MAJOR ORDER MUST BE SET HERE   #
#############################################
COLUMN_MAJOR=false
ifeq ($(COLUMN_MAJOR), true)
	CFLAGS += -DCOLUMN_MAJOR_ORDER
endif

# Should be a cleaner way to do this if we rename the exec sources
SRC=$(wildcard src/*.c)
OBJ=$(patsubst %.c,%.o,$(SRC))

.PHONY: all clean doc test cover

all: lib/libgensvm.a

ifdef NOATLAS
override LDFLAGS+=-lcblas -llapack -lm
else
override LDFLAGS+=-lcblas -llapack -lm -latlas
endif

debug: CFLAGS += -DDEBUG
debug: all

doc: cover
	$(DOXY) $(DOXYFILE)

clean:
	rm -rf *.o src/*.o lib/*.a *.{gcno,gcov} src/*.{gcno,gcda}
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

src/%.o: src/%.c
	@$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -c $< -o $@
	@echo $< ...
