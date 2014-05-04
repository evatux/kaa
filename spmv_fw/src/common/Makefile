include ../Makefile.config

COMMON          = unit core graphics
COMMON_DEPS     = $(addsuffix .h, $(COMMON)) common.h
COMMON_OBJS     = $(addprefix $(DIR.BUILD)/, $(addsuffix .o, $(COMMON)))

INC_DIRS        = $(DIR.COMMON)
INCLUDE         = $(addprefix -I, $(INC_DIRS))

all: $(COMMON_OBJS)

$(DIR.BUILD)/%.o: %.c $(COMMON_DEPS)
	$(C.COMP) -c $< -o $@ $(INCLUDE)

clean:
	-rm -f *.o
	-rm -f *~
