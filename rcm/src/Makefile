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
RCM_DIR   		= $(ROOT_DIR)/rcm/src

RCM       		= runit  rcm rmd md nd
RCM_DEPS  		= $(addsuffix .h, $(RCM))
RCM_OBJS  		= $(addsuffix .o, $(RCM))

COMMON       	= unit core graphics cholesky solver
COMMON_DEPS  	= $(addsuffix .h, $(COMMON)) common.h
COMMON_OBJS  	= $(addsuffix .o, $(COMMON))
COMMON_OBJS_WITH_PATH = $(addprefix $(COMMON_DIR)/, $(COMMON_OBJS))

all:  rcm.out

rcm.out: main.o $(RCM_OBJS)
	cd $(COMMON_DIR) && make Wall=$(Wall) DEBUG=$(DEBUG) $(COMMON_OBJS)
	$(CC) $^ $(COMMON_OBJS_WITH_PATH) $(LIBS) $(PNGLIBS) -o ../bin/$@

%.o: %.c $(RCM_DEPS)
	$(CC) -c $< $(CFLAGS) -I$(RCM_DIR) -I$(COMMON_DIR) -o $@

full_clean: clean
	cd $(COMMON_DIR) && make clean

clean:
	-rm -f *.o
	-rm -f *~
