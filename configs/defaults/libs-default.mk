override LIBS += -L$(LIBS_DIR)/2decomp-fft -ldecomp2d
override INCS += -I$(LIBS_DIR)/2decomp-fft/mod

ifeq ($(strip $(GPU)),1)
override LIBS += -L$(LIBS_DIR)/cuDecomp/build/lib -lcudecomp -lcudecomp_fort -cudalib=cufft
override INCS += -I$(LIBS_DIR)/cuDecomp/build/include
endif

ifeq ($(strip $(USE_NVTX)),1)
NVHPC_HOME ?= /opt/nvidia/hpc_sdk/Linux_x86_64/2022
override LIBS += -L$(NVHPC_HOME)/cuda/lib64 -lnvToolsExt
endif

ifneq ($(strip $(GPU)),1)
override LIBS += -lfftw3

ifeq ($(strip $(OPENMP)),1)
override LIBS += -lfftw3_threads
endif

ifeq ($(strip $(SINGLE_PRECISION)),1)
override LIBS += -lfftw3f
ifeq ($(strip $(OPENMP)),1)
override LIBS += -lfftw3f_threads
endif
endif

endif
