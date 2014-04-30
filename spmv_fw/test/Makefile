SRC_DIR = ../src
CUR_DIR = $(shell pwd)

.PHONY: $(SRC_DIR)

$(SRC_DIR):
	$(MAKE) -C $@ OUT_DIR=$(CUR_DIR) CC=icc DEFS=-O3

clean:
	$(MAKE) -C $(SRC_DIR) clean
