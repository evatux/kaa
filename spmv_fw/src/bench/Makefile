include ../Makefile.config

PNGLIB = -lpng
LIBS   = $(MKLNS.LIBS) $(MKL.LIBS) -lm -lpthread $(PNGLIB)

BENCH           =  correctness perf load_config bench
BENCH_OUT       =  bench.out
BENCH_DEPS      =  bench.h
BENCH_OBJS      =  $(addprefix $(DIR.BUILD)/, $(addsuffix .o, $(BENCH)))

ALL_OBJS        =  $(wildcard  $(DIR.BUILD)/*.o)
ALL_OBJS        += $(BENCH_OBJS)

INC_DIRS        =  $(DIR.COMMON) $(DIR.SPMV) $(DIR.KERS)
INCLUDE         =  $(addprefix -I, $(INC_DIRS))

.SECONDARY:

all:  $(DIR.OUT)/$(BENCH_OUT)

$(DIR.BUILD)/%.o: %.c $(BENCH_DEPS)
	$(C.COMP) -c $< -o $@ $(INCLUDE)

$(DIR.OUT)/%.out: $(BENCH_OBJS)
	$(C.COMP.OMP) $(ALL_OBJS) -o $@ $(LIBS)

clean:
	-rm -f $(BENCH_OBJS)
	-rm -f $(DIR.OUT)/bench.out
