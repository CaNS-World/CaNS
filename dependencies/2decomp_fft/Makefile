2DECOMP_DIR=$(CURDIR)

.PHONY: lib examples clean install_dir

all: lib basic_test

lib:
	cd lib; $(MAKE) $@

examples:
	cd $@ ; $(MAKE) $@

basic_test: examples
	@echo "Basic Test target is examples"

clean:
	cd src; $(MAKE) $@
	cd lib; $(MAKE) $@
	cd include; rm -f *.mod
	cd examples; $(MAKE) $@

install_dir:
	mkdir -p $(DESTDIR)$(prefix)
	mkdir -p $(DESTDIR)$(prefix)/include
	mkdir -p $(DESTDIR)$(prefix)/lib
	mkdir -p $(DESTDIR)$(prefix)/doc

install: all install_dir
	cp $(2DECOMP_DIR)/include/*.mod $(DESTDIR)$(prefix)/include
	cp $(2DECOMP_DIR)/lib/lib*.a $(DESTDIR)$(prefix)/lib
	cp $(2DECOMP_DIR)/README $(DESTDIR)$(prefix)/README_2DECOMP
	cp $(2DECOMP_DIR)/doc/* $(DESTDIR)$(prefix)/doc
