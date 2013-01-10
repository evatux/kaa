CC = gcc
DEPS = unit.h rcm.h
CFLAGS += -std=c99 -I.
#LIBS = -fopenmp

ifeq ($(DEBUG),1)
CFLAGS += -D_DEBUG
endif

all: a.out clean

a.out: main.o rcm.o unit.o
	$(CC) $^ $(LIBS) -o $@

%.o: %.c $(DEPS)
	$(CC) -c $< $(CFLAGS) $(LIBS) -o $@

clean:
	-rm -f *.o
	-rm -f *~
