ifeq ($(strip $(FCOMP)),GNU)
FFLAGS_MOD_DIR := -J # extra space
endif
ifeq ($(findstring INTEL,$(strip $(FCOMP))),INTEL)
FFLAGS_MOD_DIR := -module # extra space
endif
ifeq ($(strip $(FCOMP)),NVIDIA)
FFLAGS_MOD_DIR := -module # extra space
override FFLAGS += -cuda -gpu=cc70,cc80,cc90
ifeq ($(strip $(GPU)),1)
ifneq ($(filter $(GPU_BACKEND),OACC OMP),$(GPU_BACKEND))
override FFLAGS += -acc -Minfo=accel
endif
ifeq ($(strip $(GPU_BACKEND)),OACC)
override FFLAGS += -acc -Minfo=accel
endif
ifeq ($(strip $(GPU_BACKEND)),OMP)
override FFLAGS += -mp=gpu -Minfo=mp
endif
endif
endif
ifeq ($(strip $(FCOMP)),CRAY)
FFLAGS_MOD_DIR := -I./build -ef -J
ifeq ($(strip $(GPU)),1)
ifneq ($(filter $(GPU_BACKEND),OACC OMP),$(GPU_BACKEND))
override FFLAGS += -hacc -hnoomp
endif
ifeq ($(strip $(GPU_BACKEND)),OACC)
override FFLAGS += -hacc -hnoomp
endif
ifeq ($(strip $(GPU_BACKEND)),OMP)
override FFLAGS += -homp -hnoacc
endif
else
override FFLAGS += -hnoacc -hnoomp
endif
endif
ifeq ($(strip $(FCOMP)),FUJITSU)
FFLAGS_MOD_DIR := -M # extra space
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
ifeq ($(strip $(FCOMP)),FUJITSU)
override FFLAGS += -g -g0
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
ifeq ($(strip $(FCOMP)),FUJITSU)
override FFLAGS += -g -g0
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
ifeq ($(strip $(FCOMP)),FUJITSU)
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
ifeq ($(strip $(FCOMP)),FUJITSU)
override FFLAGS += -O3 -Kfast
endif
  
endif

CUSTOM_DEFINES =  SINGLE_PRECISION \
		          LOOP_UNSWITCHING \
		          DECOMP_X_IO \
				  USE_DIEZDECOMP \
		          USE_HIP \
                  USE_NVTX

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

ifeq ($(strip $(OPENMP)),1)
ifeq      ($(strip $(FCOMP)),GNU)
override FFLAGS += -fopenmp
else ifeq ($(findstring INTEL,$(strip $(FCOMP))),INTEL)
override FFLAGS += -qopenmp
else ifeq ($(strip $(FCOMP)),NVIDIA)
override FFLAGS += -mp
else ifeq ($(strip $(FCOMP)),CRAY)
override FFLAGS += -fopenmp
else ifeq ($(strip $(FCOMP)),FUJITSU)
override FFLAGS += -fopenmp
else
override FFLAGS += -fopenmp
endif
else
ifeq ($(strip $(FCOMP)),CRAY)
override FFLAGS += -fno-openmp
endif
endif
