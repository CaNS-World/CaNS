#=======================================================================
# Makefile for 2DECOMP&FFT in a CaNS checkout
#=======================================================================
# Choose pre-processing options
#   -DDOUBLE_PREC - use double precision
#   -DOVERWRITE   - enable overwriting optimisation in MKL FFT

CFG_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
CANS_DIR := $(abspath $(CFG_DIR)/../../..)

GIT_VERSION := $(shell git describe --tag --long --always)

FFT ?= generic
PARAMOD ?= mpi
PROFILER ?= none
BUILD_CONFIG_FILE ?= build.conf
BUILD ?=
FCFLAGS ?=
LDFLAGS ?=
FFLAGS = $(FCFLAGS)
LFLAGS = $(LDFLAGS)
MODFLAG = -J

LIBDECOMP = decomp2d

AR = ar
LIBOPT = rcs

CANS_CONF = $(CANS_DIR)/$(BUILD_CONFIG_FILE)
ifeq ("$(wildcard $(CANS_CONF))","")
$(error Expected CaNS build configuration at $(CANS_CONF))
endif
include $(CFG_DIR)/Makefile-CaNS.compilers

SRCDECOMP = decomp_2d_constants.f90 factor.f90 decomp_2d_mpi.f90 decomp_2d.f90 log.f90 io.f90

ifeq ($(PARAMOD),gpu)
  SRCDECOMP += decomp_2d_cumpi.f90
  ifeq ($(NCCL),yes)
    SRCDECOMP += decomp_2d_nccl.f90
  endif
endif

ifeq ($(FFT),fftw3)
  FFTW3_PATH ?= /usr
  FFTW3_PATH_INCLUDE ?= $(FFTW3_PATH)/include
  FFTW3_PATH_LIB ?= $(FFTW3_PATH)/lib
  INC = -I$(FFTW3_PATH_INCLUDE)
  LIBFFT = -L$(FFTW3_PATH_LIB) -lfftw3 -lfftw3f
else ifeq ($(FFT),fftw3_f03)
  FFTW3_PATH ?= /usr
  FFTW3_PATH_INCLUDE ?= $(FFTW3_PATH)/include
  FFTW3_PATH_LIB ?= $(FFTW3_PATH)/lib
  INC = -I$(FFTW3_PATH_INCLUDE)
  LIBFFT = -L$(FFTW3_PATH_LIB) -lfftw3 -lfftw3f
else ifeq ($(FFT),generic)
  SRCDECOMP += glassman.f90
  INC =
  LIBFFT =
else ifeq ($(FFT),mkl)
  SRCDECOMP += $(MKLROOT)/include/mkl_dfti.f90
  LIBFFT = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread
  INC = -I$(MKLROOT)/include
else ifeq ($(FFT),cufft)
  CUFFT_PATH ?= $(NVHPC)/Linux_x86_64/$(EBVERSIONNVHPC)/compilers
  INC = -I$(CUFFT_PATH)/include
endif

LIBIO :=
OPTIO :=
INCIO :=
ADIOS2DIR :=
ifeq ($(IO),adios2)
  ifeq ($(ADIOS2DIR),)
    $(error Set ADIOS2DIR=/path/to/adios2/install/)
  endif
  OPTIO := -DADIOS2 $(OPT)
  INCIO := $(INC) $(shell $(ADIOS2DIR)/bin/adios2-config --fortran-flags)
  LIBIO := $(shell $(ADIOS2DIR)/bin/adios2-config --fortran-libs)
endif

ifneq ($(PROFILER),none)
  DEFS += -DPROFILER
endif
ifeq ($(PROFILER),caliper)
  CALIPER_PATH ?= xxxxxxxxx/caliper/caliper_2.8.0
  SRCDECOMP += profiler_caliper.f90
  INC += -I$(CALIPER_PATH)/include/caliper/fortran
  LFLAGS += -L$(CALIPER_PATH)/lib -lcaliper
endif

OPT =
LINKOPT = $(FFLAGS)

OBJDIR = obj
SRCDIR = src
DECOMPINC = mod
FFLAGS += $(MODFLAG)$(DECOMPINC) -I$(DECOMPINC)

SRCDECOMP += fft_$(FFT).f90 fft_log.f90
SRCDECOMP_ = $(patsubst %.f90,$(SRCDIR)/%.f90,$(filter-out %/mkl_dfti.f90,$(SRCDECOMP)))
SRCDECOMP_ += $(filter %/mkl_dfti.f90,$(SRCDECOMP))
OBJDECOMP_MKL_ = $(patsubst $(MKLROOT)/include/%.f90,$(OBJDIR)/%.f90,$(filter %/mkl_dfti.f90,$(SRCDECOMP_)))
OBJDECOMP_MKL = $(OBJDECOMP_MKL_:%.f90=%.o)
OBJDECOMP = $(SRCDECOMP_:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)

OPT += $(OPTIO)
INC += $(INCIO)

-include Makefile.settings

all: $(DECOMPINC) $(OBJDIR) $(LIBDECOMP)

$(DECOMPINC):
	mkdir -p $(DECOMPINC)

$(LIBDECOMP): Makefile.settings lib$(LIBDECOMP).a

lib$(LIBDECOMP).a: $(OBJDECOMP_MKL) $(OBJDECOMP)
	$(AR) $(LIBOPT) $@ $^

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDECOMP): $(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@

$(OBJDECOMP_MKL): $(OBJDIR)/%.o : $(MKLROOT)/include/%.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $(MKLROOT)/include/mkl_dfti.f90 -o $(OBJDIR)/mkl_dfti.o

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o $(DECOMPINC)/*.mod $(DECOMPINC)/*.smod lib$(LIBDECOMP).a
	rm -f ./*.o ./*.mod ./*.smod
	rm -f Makefile.settings

.PHONY: Makefile.settings

Makefile.settings:
	echo "FC = $(FC)" > $@
	echo "FFLAGS = $(FFLAGS)" >> $@
	echo "OPT = $(OPT)" >> $@
	echo "DEFS = $(DEFS)" >> $@
	echo "INC = $(INC)" >> $@
	echo "LIBOPT = $(LIBOPT)" >> $@
	echo "LIBFFT = $(LIBFFT)" >> $@
	echo "LFLAGS = $(LFLAGS)" >> $@

export
