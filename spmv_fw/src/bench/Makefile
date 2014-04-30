CC = gcc
DEFS =
CFLAGS = -g -std=c99
PNGLIB = -lpng
LIBS = -lm -lpthread $(PNGLIB)

ifeq ($(DEBUG),1)
CFLAGS += -D_DEBUG
endif

ifneq ($(MKL),0)
ifndef MKLROOT
MKLROOT  = /opt/intel/composer_xe_2013_sp1.2.144/mkl
endif
MKLLIBS  = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
MKLLDIR  = $(MKLROOT)/lib/intel64 $(MKLROOT)/../compiler/lib/intel64
MKL_LIBS = $(addprefix -L, $(MKLLDIR)) -Wl,-Bstatic $(MKLLIBS) -Wl,-Bdynamic
CFLAGS  += -D_MKL_ENABLE
endif

ifneq ($(MKL_NEW_SPMV),0)
ifndef MKLNSROOT
MKLNSROOT   = /opt/intel/composer_xe_2013_sp1.2.144/mkl
endif
MKLNSLIBS  = -lmkl_spmv_pkg_sse
MKLNSLDIR  = $(MKLNSROOT)/lib/intel64
MKLNS_LIBS = $(addprefix -L, $(MKLNSLDIR)) -Wl,-Bstatic $(MKLNSLIBS) -Wl,-Bdynamic
CFLAGS     += -D_MKL_NEW_SPMV_ENABLE
endif

ifeq ($(Wall),1)
CFLAGS += -Wall
endif

ROOT_DIR        = ..
COMMON_DIR      = $(ROOT_DIR)/common
SPMV_DIR        = $(ROOT_DIR)/spmv
KERS_DIR        = $(SPMV_DIR)/kernels
BENCH_DIR       = $(ROOT_DIR)/bench
OUT_DIR         = .

BENCH           = correctness perf load_config bench
BENCH_OUT       = bench.out
BENCH_DEPS      = bench.h
BENCH_OBJS      = $(addsuffix .o, $(BENCH))

ALL_DIRS        = $(COMMON_DIR) $(SPMV_DIR) $(KERS_DIR)
ALL_OBJS        = $(foreach dir, $(ALL_DIRS), $(wildcard $(dir)/*.o))
ALL_OBJS        += $(BENCH_OBJS)

INC_DIR         += $(COMMON_DIR) $(SPMV_DIR) $(KERS_DIR)

all:  $(BENCH_OUT)

%.o: %.c $(BENCH_DEPS)
	$(CC) $(CFLAGS) $(DEFS) -c $< $(addprefix -I, $(INC_DIR)) -o $@

%.out: $(BENCH_OBJS)
	$(CC) $(CFLAGS) $(DEFS) $(ALL_OBJS) -o $(OUT_DIR)/$@ $(MKLNS_LIBS) $(MKL_LIBS) $(LIBS)

clean:
	-rm -f *.out
	-rm -f *.o
	-rm -f *~
