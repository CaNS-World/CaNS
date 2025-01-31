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

CUSTOM_DEFINES =  DEBUG \
				  DEBUG_SOLVER \
				  TIMING \
				  IMPDIFF IMPDIFF_1D \
				  DECOMP_X DECOMP_Y DECOMP_Z \
				  SINGLE_PRECISION \
				  DECOMP_X_IO \
                  USE_NVTX \
				  GRIDPOINT_NATURAL_CHANNEL \
				  MASK_DIVERGENCE_CHECK \
				  BOUSSINESQ_BUOYANCY

DEFINES += $(foreach var,$(CUSTOM_DEFINES),$(if $(filter 1,$(strip $($(var)))), -D_$(var)))

ifeq ($(strip $(IMPDIFF_1D)),1)
DEFINES += -D_IMPDIFF -D_IMPDIFF_1D
endif
ifeq      ($(strip $(PENCIL_AXIS)),1)
DEFINES += -D_DECOMP_X
else ifeq ($(strip $(PENCIL_AXIS)),2)
DEFINES += -D_DECOMP_Y
else ifeq ($(strip $(PENCIL_AXIS)),3)
DEFINES += -D_DECOMP_Z
endif

DEFINES := $(sort $(DEFINES)) # remove duplicates

ifeq ($(strip $(BOUSSINESQ_BUOYANCY)),1)
DEFINES += -D_BOUSSINESQ_BUOYANCY
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
