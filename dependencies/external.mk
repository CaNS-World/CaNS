#
# external libraries compilation
#
D2D_MAKEFILE := ../configs/2decomp-fft/Makefile
DIEZ_MAKEFILE := ../configs/diezDecomp/Makefile

ifeq ($(strip $(GPU)),1)
ifeq ($(strip $(USE_DIEZDECOMP)),1)
libs: $(wildcard $(LIBS_DIR)/2decomp-fft/src/*.f90)
	$(MAKE) -C $(LIBS_DIR)/2decomp-fft -f $(D2D_MAKEFILE)
	$(MAKE) -C $(LIBS_DIR)/diezDecomp -f $(DIEZ_MAKEFILE) -j
libsclean: $(wildcard $(LIBS_DIR)/2decomp-fft/src/*.f90)
	$(MAKE) -C $(LIBS_DIR)/2decomp-fft -f $(D2D_MAKEFILE) clean
	$(MAKE) -C $(LIBS_DIR)/diezDecomp -f $(DIEZ_MAKEFILE) clean
else
ifeq ($(strip $(USE_HIP)),1)
HIP_GPU_TARGET ?= native
libs: $(wildcard $(LIBS_DIR)/2decomp-fft/src/*.f90)
	$(MAKE) -C $(LIBS_DIR)/2decomp-fft -f $(D2D_MAKEFILE)
	cd $(LIBS_DIR)/hipDecomp && mkdir -p build && cd build && cmake .. -DHIPDECOMP_USE_CUDECOMP_NAMES=1 -DCMAKE_HIP_ARCHITECTURES="$(HIP_GPU_TARGET)" && make -j
libsclean: $(wildcard $(LIBS_DIR)/2decomp-fft/src/*.f90)
	$(MAKE) -C $(LIBS_DIR)/2decomp-fft -f $(D2D_MAKEFILE) clean
	cd $(LIBS_DIR)/hipDecomp/build && make clean; cd .. && rm -rf build
else
libs: $(wildcard $(LIBS_DIR)/2decomp-fft/src/*.f90)
	$(MAKE) -C $(LIBS_DIR)/2decomp-fft -f $(D2D_MAKEFILE)
	cd $(LIBS_DIR)/cuDecomp && mkdir -p build && cd build && cmake .. && make -j
libsclean: $(wildcard $(LIBS_DIR)/2decomp-fft/src/*.f90)
	$(MAKE) -C $(LIBS_DIR)/2decomp-fft -f $(D2D_MAKEFILE) clean
	cd $(LIBS_DIR)/cuDecomp/build && make clean; cd .. && rm -rf build
endif
endif
else
libs: $(wildcard $(LIBS_DIR)/2decomp-fft/src/*.f90)
	$(MAKE) -C $(LIBS_DIR)/2decomp-fft -f $(D2D_MAKEFILE)
libsclean: $(wildcard $(LIBS_DIR)/2decomp-fft/src/*.f90)
	$(MAKE) -C $(LIBS_DIR)/2decomp-fft -f $(D2D_MAKEFILE) clean
endif
