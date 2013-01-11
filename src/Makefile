CC = gcc
DEPS = unit.h rcm.h
CFLAGS += -std=c99 -I.
#LIBS = -fopenmp
PNGLIBS = -lpng

ifeq ($(DEBUG),1)
CFLAGS += -D_DEBUG
endif

all: rcm.out

rcm.out: main.o rcm.o core.o graphics.o unit.o
	$(CC) $^ $(LIBS) $(PNGLIBS) -o ../bin/$@

%.o: %.c $(DEPS)
	$(CC) -c $< $(CFLAGS) $(LIBS) -o $@

clean:
	-rm -f *.o
	-rm -f *~
