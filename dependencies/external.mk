library: $(wildcard $(LIBS_DIR)/2decomp_fft/src/*.f90)
	cd $(LIBS_DIR)/2decomp_fft && make clean && make
libclean: $(wildcard $(LIBS_DIR)/2decomp_fft/src/*.f90)
	cd $(LIBS_DIR)/2decomp_fft && make clean
