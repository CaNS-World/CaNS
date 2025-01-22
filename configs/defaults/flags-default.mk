ifeq ($(strip $(FCOMP)),GNU)
FFLAGS_MOD_DIR := -J
endif
ifeq ($(findstring INTEL,$(strip $(FCOMP))),INTEL)
FFLAGS_MOD_DIR := -module
endif
ifeq ($(strip $(FCOMP)),NVIDIA)
FFLAGS_MOD_DIR := -module
ifeq ($(strip $(GPU)),1)
override FFLAGS += -acc -cuda -Minfo=accel -gpu=cc60,cc70,cc80
endif
endif
ifeq ($(strip $(FCOMP)),CRAY)
FFLAGS_MOD_DIR := -J
endif

ifeq ($(strip $(FFLAGS_DEBUG)),1)

ifeq ($(strip $(FCOMP)),GNU)
override FFLAGS += -O0 -g -fbacktrace -Wall -Wextra -pedantic -fcheck=all -finit-real=snan -ffpe-trap=invalid -std=f2018
endif
ifeq ($(findstring INTEL,$(strip $(FCOMP))),INTEL)
override FFLAGS += -O0 -g -traceback -fpe0 -stand f18
endif
ifeq ($(strip $(FCOMP)),NVIDIA)
override FFLAGS += -O0 -g -traceback -Mstandard -Minform=inform -Mbackslash -Mbounds -Mchkstk
endif
ifeq ($(strip $(FCOMP)),CRAY)
override FFLAGS += -g -G0
endif
  
endif

ifeq ($(strip $(FFLAGS_DEBUG_MAX)),1)

ifeq ($(strip $(FCOMP)),GNU)
override FFLAGS += -O0 -g -fbacktrace -Wall -Wextra -Wimplicit-interface -Wno-unused-function -fPIC -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-99999999 -std=f2018
endif
ifeq ($(findstring INTEL,$(strip $(FCOMP))),INTEL)
override FFLAGS += -O0 -warn all -g -traceback -fpe0 -stand f18
endif
ifeq ($(strip $(FCOMP)),NVIDIA)
override FFLAGS += -O0 -g -traceback -Ktrap=fp -Mstandard -Minform=inform -Mbackslash -Mbounds -Mchkptr -Mchkstk
endif
ifeq ($(strip $(FCOMP)),CRAY)
override FFLAGS += -g -G0
endif

endif

ifeq ($(strip $(FFLAGS_OPT)),1)

ifeq ($(strip $(FCOMP)),GNU)
override FFLAGS += -O3
endif
ifeq ($(findstring INTEL,$(strip $(FCOMP))),INTEL)
override FFLAGS += -O3
endif
ifeq ($(strip $(FCOMP)),NVIDIA)
override FFLAGS += -O3 -fast
endif
ifeq ($(strip $(FCOMP)),CRAY)
override FFLAGS += -O3
endif
  
endif

ifeq ($(strip $(FFLAGS_OPT_MAX)),1)

ifeq ($(strip $(FCOMP)),GNU)
override FFLAGS += -Ofast -march=native
endif
ifeq ($(findstring INTEL,$(strip $(FCOMP))),INTEL)
override FFLAGS += -fast -xHost
endif
ifeq ($(strip $(FCOMP)),NVIDIA)
override FFLAGS += -O3 -fast -Mnouniform -Mfprelaxed
endif
ifeq ($(strip $(FCOMP)),CRAY)
override FFLAGS += -O3 -hfp3
endif
  
endif

ifeq ($(strip $(DEBUG)),1)
DEFINES += -D_DEBUG
endif
ifeq ($(strip $(DEBUG_SOLVER)),1)
DEFINES += -D_DEBUG_SOLVER
endif
ifeq ($(strip $(TIMING)),1)
DEFINES += -D_TIMING
endif
ifeq ($(strip $(IMPDIFF)),1)
DEFINES += -D_IMPDIFF
endif
ifeq ($(strip $(IMPDIFF_1D)),1)
DEFINES += -D_IMPDIFF -D_IMPDIFF_1D
endif
ifeq      ($(strip $(DECOMP_X)),1)
DEFINES += -D_DECOMP_X
else ifeq ($(strip $(DECOMP_Y)),1)
DEFINES += -D_DECOMP_Y
else ifeq ($(strip $(DECOMP_Z)),1)
DEFINES += -D_DECOMP_Z
endif
ifeq      ($(strip $(PENCIL_AXIS)),1)
DEFINES += -D_DECOMP_X
else ifeq ($(strip $(PENCIL_AXIS)),2)
DEFINES += -D_DECOMP_Y
else ifeq ($(strip $(PENCIL_AXIS)),3)
DEFINES += -D_DECOMP_Z
endif
ifeq ($(strip $(SINGLE_PRECISION)),1)
DEFINES += -D_SINGLE_PRECISION
endif

ifeq      ($(strip $(DECOMP_X_IO)),1)
DEFINES += -D_DECOMP_X_IO
endif

ifeq ($(strip $(USE_NVTX)),1)
DEFINES += -D_USE_NVTX
endif

ifeq ($(strip $(GRIDPOINT_NATURAL_CHANNEL)),1)
DEFINES += -D_GRIDPOINT_NATURAL_CHANNEL
endif
ifeq ($(strip $(MASK_DIVERGENCE_CHECK)),1)
DEFINES += -D_MASK_DIVERGENCE_CHECK
endif

ifeq ($(strip $(OPENMP)),1)
ifeq      ($(strip $(FCOMP)),GNU)
override FFLAGS += -fopenmp
else ifeq ($(findstring INTEL,$(strip $(FCOMP))),INTEL)
override FFLAGS += -qopenmp
else ifeq ($(strip $(FCOMP)),NVIDIA)
override FFLAGS += -mp
else ifeq ($(strip $(FCOMP)),CRAY)
override FFLAGS += -homp
else
override FFLAGS += -fopenmp
endif
endif
