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
FEM_DIR			= $(ROOT_DIR)/fem

FEM				= 
FEM_DEPS		= $(addsuffix .h, $(RCM))
FEM_OBJS		= $(addsuffix .o, $(RCM))

COMMON			= unit core graphics
COMMON_DEPS		= $(addsuffix .h, $(COMMON)) common.h
COMMON_OBJS		= $(addsuffix .o, $(COMMON))
COMMON_OBJS_WITH_PATH = $(addprefix $(COMMON_DIR)/, $(COMMON_OBJS))

all:  laplace_2d.out laplace_3d.out

common:
	cd $(COMMON_DIR) && make Wall=$(Wall) DEBUG=$(DEBUG) $(COMMON_OBJS)

%.out: %.o common
	$(CC) $< $(COMMON_OBJS_WITH_PATH) $(LIBS) $(PNGLIBS) -o $@

%.o: %.c $(FEM_DEPS)
	$(CC) -c $< $(CFLAGS) -I$(FEM_DIR) -I$(COMMON_DIR) -o $@

full_clean: clean
	cd $(COMMON_DIR) && make clean

clean:
	-rm -f *.o
	-rm -f *~
