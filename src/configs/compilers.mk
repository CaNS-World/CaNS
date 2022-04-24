FC = mpifort
ifeq ($(strip $(FCOMP)),GNU)
FC = mpifort
endif
ifeq ($(strip $(FCOMP)),INTEL)
FC = mpiifort
endif
ifeq ($(strip $(FCOMP)),NVIDIA)
FC = mpifort
endif
