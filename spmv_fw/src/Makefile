include Makefile.config

SUBDIRS = common spmv bench

.PHONY: subdirs $(SUBDIRS)

subdirs: prepare $(SUBDIRS)

prepare:
	-mkdir -p $(DIR.BUILD)

$(SUBDIRS):
	$(MAKE) -C $@

common: prepare
spmv:   common
bench:  spmv

clean:
	-rm -f $(DIR.BUILD)/*.o
	-rm -f $(DIR.OUT)/bench.out
