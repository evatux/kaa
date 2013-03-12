CC = gcc
DEPS = unit.h core.h rcm.h md.h nd.h graphics.h cholesky.h solver.h
CFLAGS += -g -std=c99 -I.
LIBS = -lm
PNGLIBS = -lpng

ifeq ($(DEBUG),1)
CFLAGS += -D_DEBUG
endif

ifeq ($(Wall),1)
CFLAGS += -Wall
endif

all:  rcm.out

test.out: test.o core.o rcm.o md.o nd.o graphics.o unit.o cholesky.o solver.o
	$(CC) $^ $(LIBS) $(PNGLIBS) -o ../bin/$@

rcm.out: main.o core.o rcm.o md.o nd.o graphics.o unit.o cholesky.o solver.o
	$(CC) $^ $(LIBS) $(PNGLIBS) -o ../bin/$@

%.o: %.c $(DEPS)
	$(CC) -c $< $(CFLAGS) -o $@

clean:
	-rm -f *.o
	-rm -f *~
