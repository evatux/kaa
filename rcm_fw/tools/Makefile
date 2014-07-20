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

ROOT_DIR		= ..
COMMON_DIR		= $(ROOT_DIR)/common
TOOLS_DIR		= $(ROOT_DIR)/tools

TOOLS			=
TOOLS_DEPS		= $(addsuffix .h, $(TOOLS))
TOOLS_OBJS		= $(addsuffix .o, $(TOOLS))
TOOLS_LIBS_DIR  = $(TOOLS_DIR)/libs
TOOLS_INC_DIR   = $(TOOLS_DIR)/include

COMMON			= unit core graphics
COMMON_DEPS		= $(addsuffix .h, $(COMMON)) common.h
COMMON_OBJS		= $(addsuffix .o, $(COMMON))
COMMON_OBJS_WITH_PATH = $(addprefix $(COMMON_DIR)/, $(COMMON_OBJS))

all:  reorder_metis.out

common:
	cd $(COMMON_DIR) && make Wall=$(Wall) DEBUG=$(DEBUG) $(COMMON_OBJS)

reorder_metis.out: reorder_metis.o common
	$(CC) $< $(COMMON_OBJS_WITH_PATH) -L$(TOOLS_LIBS_DIR) -lmetis $(LIBS) $(PNGLIBS) -o $@

%.o: %.c $(TOOLS_DEPS)
	$(CC) -c $< $(CFLAGS) -I$(TOOLS_DIR) -I$(TOOLS_INC_DIR) -I$(COMMON_DIR) -o $@

full_clean: clean
	cd $(COMMON_DIR) && make clean

clean:
	-rm -f *.o
	-rm -f *~
