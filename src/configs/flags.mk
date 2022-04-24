ifeq ($(strip $(FCOMP)),GNU)
FFLAGS_MOD_DIR := -J
endif
ifeq ($(strip $(FCOMP)),INTEL)
FFLAGS_MOD_DIR := -module
endif
ifeq ($(strip $(FCOMP)),NVIDIA)
FFLAGS_MOD_DIR := -module
endif

ifeq ($(strip $(FFLAGS_DEBUG)),1)

ifeq ($(strip $(FCOMP)),GNU)
FFLAGS += -O0 -g -fbacktrace -Wall -Wextra -pedantic -fcheck=all -finit-real=snan -ffpe-trap=invalid -std=f2018
endif

ifeq ($(strip $(FCOMP)),INTEL)
FFLAGS += -O0 -g -traceback -fpe0 -stand f18
endif

ifeq ($(strip $(FCOMP)),NVIDIA)
FFLAGS += -O0 -g -traceback -Mstandard -Minform=inform -Mbackslash -Mbounds -Mchkptr -Mchkstk
endif
  
endif
ifeq ($(strip $(FFLAGS_DEBUG_MAX)),1)

ifeq ($(strip $(FCOMP)),GNU)
FFLAGS += -O0 -Wall -Wextra -Wimplicit-interface -Wno-unused-function -fPIC -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-99999999 -std=f2018
endif

endif

ifeq ($(strip $(FFLAGS_OPT)),1)

ifeq ($(strip $(FCOMP)),GNU)
FFLAGS += -O3
endif

ifeq ($(strip $(FCOMP)),INTEL)
FFLAGS += -O3
endif

ifeq ($(strip $(FCOMP)),NVIDIA)
FFLAGS += -O3
endif
  
endif

ifeq ($(strip $(FFLAGS_OPT_MAX)),1)

ifeq ($(strip $(FCOMP)),GNU)
FFLAGS += -Ofast -march=native
endif

ifeq ($(strip $(FCOMP)),INTEL)
FFLAGS += -fast -xHost
endif

ifeq ($(strip $(FCOMP)),NVIDIA)
FFLAGS += -fast
endif
  
endif

ifeq ($(strip $(DEBUG)),1)
DEFINES += -D_DEBUG
endif
ifeq ($(strip $(IMPDIFF)),1)
DEFINES += -D_IMPDIFF
endif
ifeq ($(strip $(IMPDIFF_1D)),1)
DEFINES += -D_IMPDIFF_1D
endif
ifeq ($(strip $(IMPDIFF_X)),1)
DEFINES += -D_DECOMP_X
endif
ifeq ($(strip $(IMPDIFF_Y)),1)
DEFINES += -D_DECOMP_Y
endif
ifeq ($(strip $(IMPDIFF_Z)),1)
DEFINES += -D_DECOMP_Z
endif
ifeq ($(strip $(SINGLE_PRECISION)),1)
DEFINES += -D_SINGLE_PRECISION
endif
