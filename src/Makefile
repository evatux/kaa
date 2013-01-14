CC = gcc
DEPS = unit.h rcm.h
CFLAGS += -std=c99 -I.
LIBS = -lm
PNGLIBS = -lpng

ifeq ($(DEBUG),1)
CFLAGS += -D_DEBUG
endif

ifeq ($(Wall),1)
CFLAGS += -Wall
endif

all: rcm.out

rcm.out: main.o rcm.o core.o graphics.o unit.o cholesky.o
	$(CC) $^ $(LIBS) $(PNGLIBS) -o ../bin/$@

%.o: %.c $(DEPS)
	$(CC) -c $< $(CFLAGS) -o $@

clean:
	-rm -f *.o
	-rm -f *~
