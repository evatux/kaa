include ../../Makefile.config

KERS            = ref ref_omp gblock_omp $(MKL.KERS) $(MKLNS.KERS)
KERS_DEPS       = $(addsuffix .h, $(KERS))
KERS_OBJS       = $(addprefix $(DIR.BUILD)/, $(addsuffix .o, $(KERS)))

INC_DIRS        = $(DIR.COMMON) $(DIR.SPMV) $(DIR.KERS)
INCLUDE         = $(addprefix -I, $(INC_DIRS) $(MKL.IDIR) $(MKLNS.IDIR))

all:  $(KERS_OBJS)

$(DIR.BUILD)/%_omp.o: %_omp.c $(KERS_DEPS)
	$(C.COMP.OMP) -c $< -o $@ $(INCLUDE)

$(DIR.BUILD)/%.o: %.c $(KERS_DEPS)
	$(C.COMP)     -c $< -o $@ $(INCLUDE)

clean:
	-rm -f $(KER_OBJS)
