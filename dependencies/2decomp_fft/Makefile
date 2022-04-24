2DECOMP_DIR=$(CURDIR)

.PHONY: lib examples clean install_dir

all: lib

lib:
	cd lib; $(MAKE) $@

clean:
	cd src; $(MAKE) $@
	cd lib; $(MAKE) $@
	cd include; rm -f *.mod

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
