SRC_DIR = ../src
CUR_DIR = $(shell pwd)

MKL_NS_DIR = /nfs/ins/proj/mkl/perfm2/efomenko/dump/nonwork/kaa/spmv/mkl_spmv_pkg_03a/__release_lnx

.PHONY: $(SRC_DIR)

$(SRC_DIR):
	$(MAKE) -C $@ DIR.OUT=$(CUR_DIR) CC=icc CFLAGS="-O3" MKLNSROOT=$(MKL_NS_DIR)

clean:
	$(MAKE) -C $(SRC_DIR) clean DIR.OUT=$(CUR_DIR)
