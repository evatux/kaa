CC = gcc
DEFS =
CFLAGS = -g -std=c99
OMPFLAGS = -fopenmp
INC_DIR =
MKL_KERS =

ifneq ($(MKL),0)
MKLROOT  = /opt/intel/composer_xe_2013_sp1.2.144/mkl
MKLIDIR  = $(MKLROOT)/include
INC_DIR += $(MKLIDIR)
MKL_KERS = mkl_spmv
CFLAGS  += -D_MKL_ENABLE
endif

ifeq ($(DEBUG),1)
CFLAGS += -D_DEBUG
endif

ifeq ($(Wall),1)
CFLAGS += -Wall
endif

ROOT_DIR        = ../..
COMMON_DIR      = $(ROOT_DIR)/common
SPMV_DIR        = $(ROOT_DIR)/spmv
KERS_DIR        = $(SPMV_DIR)/kernels

KERS            = ref ref_omp $(MKL_KERS)
KERS_DEPS       = $(addsuffix .h, $(KERS))
KERS_OBJS       = $(addsuffix .o, $(KERS))

INC_DIR         += $(COMMON_DIR) $(SPMV_DIR) $(KERS_DIR)

all:  $(KERS_OBJS)

%_omp.o: %_omp.c $(KERS_DEPS)
	$(CC) -c $< $(OMPFLAGS) $(CFLAGS) $(DEFS) $(addprefix -I, $(INC_DIR)) -o $@

%.o: %.c $(KERS_DEPS)
	$(CC) -c $< $(CFLAGS) $(DEFS) $(addprefix -I, $(INC_DIR)) -o $@

clean:
	-rm -f *.o
	-rm -f *~
