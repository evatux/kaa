SRC_DIR = ../src
CUR_DIR = $(shell pwd)

MKL_NS_DIR = /nfs/ins/proj/mkl/perfm2/efomenko/dump/nonwork/kaa/spmv/mkl_spmv_pkg_03a/__release_lnx

.PHONY: $(SRC_DIR)

$(SRC_DIR):
	$(MAKE) -C $@ OUT_DIR=$(CUR_DIR) CC=icc DEFS="-O3 -D_DOUBLE" MKLNSROOT=$(MKL_NS_DIR)

clean:
	$(MAKE) -C $(SRC_DIR) clean
