CC = gcc
DEFS =
CFLAGS = -g -std=c99

ifeq ($(DEBUG),1)
CFLAGS += -D_DEBUG
endif

ifeq ($(Wall),1)
CFLAGS += -Wall
endif

ROOT_DIR        = ..
COMMON_DIR      = $(ROOT_DIR)/common
SPMV_DIR        = $(ROOT_DIR)/spmv
KERS_DIR        = $(SPMV_DIR)/kernels

SPMV            = spmv
SPMV_DEPS       = $(addsuffix .h, $(SPMV))
SPMV_OBJS       = $(addsuffix .o, $(SPMV))

INC_DIR         += $(COMMON_DIR) $(SPMV_DIR)

.PHONY: $(KERS_DIR)

all: $(SPMV_OBJS) $(KERS_DIR)

$(KERS_DIR):
	$(MAKE) -C $@

%.o: %.c $(SPMV_DEPS)
	$(CC) -c $< $(CFLAGS) $(DEFS) $(addprefix -I, $(INC_DIR)) -o $@

clean:
	$(MAKE) -C $(KERS_DIR) clean
	-rm -f *.o
	-rm -f *~
