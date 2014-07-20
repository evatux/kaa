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
RCM_DIR         = $(ROOT_DIR)/rcm/src

TOOLS			=
TOOLS_DEPS		= $(addsuffix .h, $(TOOLS))
TOOLS_OBJS		= $(addsuffix .o, $(TOOLS))
TOOLS_LIBS_DIR  = $(TOOLS_DIR)/libs
TOOLS_INC_DIR   = $(TOOLS_DIR)/include

RCM             = tnd runit
RCM_DEPS        = $(addsuffix .h, $(RCM))
RCM_OBJS        = $(addsuffix .o, $(RCM))
RCM_OBJS_WITH_PATH = $(addprefix $(RCM_DIR)/, $(RCM_OBJS))

COMMON			= unit core graphics
COMMON_DEPS		= $(addsuffix .h, $(COMMON)) common.h
COMMON_OBJS		= $(addsuffix .o, $(COMMON))
COMMON_OBJS_WITH_PATH = $(addprefix $(COMMON_DIR)/, $(COMMON_OBJS))

INC             = $(addprefix -I, $(TOOLS_DIR) $(TOOLS_INC_DIR) $(COMMON_DIR) $(RCM_DIR))

all:  reorder_metis.out reorder_tnd.out

common:
	cd $(COMMON_DIR) && make Wall=$(Wall) DEBUG=$(DEBUG) $(COMMON_OBJS)

rcm:
	cd $(RCM_DIR) && make Wall=$(Wall) DEBUG=$(DEBUG) $(RCM_OBJS)

reorder_metis.out: reorder_metis.o common
	$(CC) $< $(COMMON_OBJS_WITH_PATH)                       -L$(TOOLS_LIBS_DIR) -lmetis $(LIBS) $(PNGLIBS) -o $@

reorder_tnd.out: reorder_tnd.o common rcm
	$(CC) $< $(COMMON_OBJS_WITH_PATH) $(RCM_OBJS_WITH_PATH) -L$(TOOLS_LIBS_DIR) $(LIBS) $(PNGLIBS) -o $@

%.o: %.c $(TOOLS_DEPS)
	$(CC) -c $< $(CFLAGS) $(INC) -o $@

full_clean: clean
	cd $(COMMON_DIR) && make clean

clean:
	-rm -f *.o
	-rm -f *~
