CC = gcc
MYFALGS = 
CFLAGS += -g -std=c99 $(MYFLAGS)
LIBS = -lm
PNGLIBS = -lpng

ifeq ($(DEBUG),1)
CFLAGS += -D_DEBUG
endif

ifeq ($(Wall),1)
CFLAGS += -Wall
endif

ROOT_DIR  		= ../..
COMMON_DIR  	= $(ROOT_DIR)/common

COMMON       	= unit core graphics cholesky solver
COMMON_DEPS  	= $(addsuffix .h, $(COMMON)) common.h
COMMON_OBJS  	= $(addsuffix .o, $(COMMON))

all: $(COMMON_OBJS)

test.out: test.o $(COMMON_OBJS)
	$(CC) $^ $(LIBS) $(PNGLIBS) -o $@

%.o: %.c $(RCM_DEPS)
	$(CC) -c $< $(CFLAGS) -I$(COMMON_DIR) -o $@

clean:
	-rm -f *.o
	-rm -f *~
