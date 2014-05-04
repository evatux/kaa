include ../Makefile.config

SPMV            = spmv
SPMV_DEPS       = $(addsuffix .h, $(SPMV))
SPMV_OBJS       = $(addprefix $(DIR.BUILD)/, $(addsuffix .o, $(SPMV)))

INC_DIRS        = $(DIR.COMMON) $(DIR.SPMV)
INCLUDE         = $(addprefix -I, $(INC_DIRS))

.PHONY: $(DIR.KERS)

all: $(SPMV_OBJS) $(DIR.KERS)

$(DIR.KERS):
	$(MAKE) -C $@

$(DIR.BUILD)/%.o: %.c $(SPMV_DEPS)
	$(C.COMP) -c $< -o $@ $(INCLUDE)

clean:
	$(MAKE) -C $(DIR.KERS) clean
	-rm -f $(SPMV_OBJS)
