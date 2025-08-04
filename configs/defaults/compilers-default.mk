FC = mpifort
ifeq ($(strip $(FCOMP)),GNU)
FC = mpifort
endif
ifeq ($(strip $(FCOMP)),INTEL)
FC = mpiifx
endif
ifeq ($(strip $(FCOMP)),INTEL_IFORT)
FC = mpiifort
endif
ifeq ($(strip $(FCOMP)),NVIDIA)
FC = mpifort
endif
ifeq ($(strip $(FCOMP)),CRAY)
ifneq ($(strip $(GPU)),1)
FC = ftn
CPP = -eZ
else
FC = hipfc -v
endif
endif
ifeq ($(strip $(FCOMP)),FUJITSU)
FC = mpifrtpx
endif
ifeq ($(strip $(FTN_MPI_WRAPPER)),1)
FC = ftn
endif
