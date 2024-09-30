# Update this CONFIGFILE for your system or use command 'make CONFIGFILE=<path to your config file>'
PWD=$(shell pwd)
CONFIGFILE=${PWD}/configs/nvhpcsdk.conf
include ${CONFIGFILE}

MPICXX ?= ${MPI_HOME}/bin/mpicxx
MPIF90 ?= ${MPI_HOME}/bin/mpif90
NVCC ?= ${CUDA_HOME}/bin/nvcc

CXXFLAGS=-O3 -std=c++14
# Note: compiling autotune.cc with -O1 due to double free issues otherwise
CXXFLAGS_O1=-O1 -std=c++14
NVFLAGS=-O3 -std=c++14 --ptxas-options=-v -cudart shared

CUDA_CC_LIST ?= 60 70 80
CUDA_CC_LAST := $(lastword $(sort ${CUDA_CC_LIST}))
NVFLAGS += $(foreach CC,${CUDA_CC_LIST},-gencode=arch=compute_${CC},code=sm_${CC})
NVFLAGS += -gencode=arch=compute_${CUDA_CC_LAST},code=compute_${CUDA_CC_LAST}

BUILDDIR=${PWD}/build
TESTDIR=${PWD}/tests
EXAMPLEDIR=${PWD}/examples
BENCHMARKDIR=${PWD}/benchmark

OBJ=${BUILDDIR}/cudecomp.o ${BUILDDIR}/cudecomp_kernels.o ${BUILDDIR}/cudecomp_kernels_rdc.o ${BUILDDIR}/autotune.o
TESTS_CC=${TESTDIR}/cc/test ${TESTDIR}/cc/fft_test
TESTS_F90=${TESTDIR}/fortran/test ${TESTDIR}/fortran/fft_test

CUDECOMPLIB=${BUILDDIR}/lib/libcudecomp.so
CUDECOMPFLIB=${BUILDDIR}/lib/libcudecomp_fort.so
CUDECOMPMOD=${BUILDDIR}/cudecomp_m.o

INCLUDES = -I${PWD}/include -I${MPI_HOME}/include -I${CUDA_HOME}/include -I${NCCL_HOME}/include  -I${CUTENSOR_HOME}/include -I${CUDACXX_HOME}/include
LIBS = -L${CUDA_HOME}/lib64 -L${CUTENSOR_HOME}/lib64 -L${NCCL_HOME}/lib -lnccl -lcutensor -lcudart
FLIBS = -cudalib=nccl,cutensor -lstdc++ -L${CUDA_HOME}/lib64
BUILD_LIBS = -L${BUILDDIR}/lib -lcudecomp -L${CUDA_HOME}/lib64 -L${CUFFT_HOME}/lib64 -lcudart -lcufft
BUILD_FLIBS = -L${BUILDDIR}/lib -lcudecomp -lcudecomp_fort -cudalib=cufft
BUILD_INCLUDES = -I${BUILDDIR}/include -I${CUDA_HOME}/include -I${CUFFT_HOME}/include -I${PWD}/include -I${NCCL_HOME}/include

ifeq ($(strip $(ENABLE_NVSHMEM)),1)
DEFINES += -DENABLE_NVSHMEM -DNVSHMEM_USE_NCCL
INCLUDES += -I${NVSHMEM_HOME}/include
LIBS += -L${CUDA_HOME}/lib64/stubs -lcuda -lnvidia-ml
FLIBS += -L${CUDA_HOME}/lib64/stubs -lcuda -lnvidia-ml
BUILD_LIBS += -L${CUDA_HOME}/lib64/stubs -lcuda -lnvidia-ml
BUILD_FLIBS += -L${CUDA_HOME}/lib64/stubs -lcuda -lnvidia-ml
ifneq ("$(wildcard ${NVSHMEM_HOME}/lib/libnvshmem_host.so)","")
LIBS += -L${NVSHMEM_HOME}/lib -lnvshmem_host
STATIC_LIBS += ${NVSHMEM_HOME}/lib/libnvshmem_device.a
else
STATIC_LIBS += ${NVSHMEM_HOME}/lib/libnvshmem.a
endif
endif
ifeq ($(strip $(ENABLE_NVTX)),1)
DEFINES += -DENABLE_NVTX
endif

ifeq ($(strip $(MPIF90)),ftn)
DEFINES += -DMPICH
endif

LIBS += ${EXTRA_LIBS}
FLIBS += ${EXTRA_LIBS}
DEFINES += ${EXTRA_DEFINES}

LIBTARGETS = ${CUDECOMPLIB}
ifeq ($(strip $(BUILD_FORTRAN)),1)
LIBTARGETS += ${CUDECOMPFLIB}
endif

export LIBS FLIBS BUILD_LIBS BUILD_FLIBS INCLUDES BUILD_INCLUDES DEFINES MPICXX MPIF90 NVCC CXXFLAGS NVFLAGS BUILD_FORTRAN

.PHONY: all lib tests benchmark examples

all: lib tests benchmark examples

lib: ${LIBTARGETS}
	@mkdir -p ${BUILDDIR}/include
	cp ./include/cudecomp.h ${BUILDDIR}/include

tests: lib
	cd ${TESTDIR}; make CONFIGFILE=${CONFIGFILE}

benchmark: lib
	cd ${BENCHMARKDIR}; make CONFIGFILE=${CONFIGFILE}

examples: lib
	cd ${EXAMPLEDIR}; make CONFIGFILE=${CONFIGFILE}

${BUILDDIR}/autotune.o: src/autotune.cc  include/*.h include/internal/*.h
	@mkdir -p ${BUILDDIR}
	${MPICXX} -fPIC ${DEFINES} ${CXXFLAGS_O1} ${INCLUDES} -c -o $@ $<

${BUILDDIR}/%.o: src/%.cc  include/*.h include/internal/*.h
	@mkdir -p ${BUILDDIR}
	${MPICXX} -fPIC ${DEFINES} ${CXXFLAGS} ${INCLUDES} -c -o $@ $<

${BUILDDIR}/cudecomp_kernels_rdc.o: src/cudecomp_kernels_rdc.cu  include/internal/*.cuh
	@mkdir -p ${BUILDDIR}
	${NVCC} -rdc=true -Xcompiler -fPIC ${DEFINES} ${NVFLAGS} ${INCLUDES} -c -o $@ $<

${BUILDDIR}/%.o: src/%.cu  include/internal/*.cuh
	@mkdir -p ${BUILDDIR}
	${NVCC} -Xcompiler -fPIC ${DEFINES} ${NVFLAGS} ${INCLUDES} -c -o $@ $<

${CUDECOMPMOD}: src/cudecomp_m.cuf 
	@mkdir -p ${BUILDDIR}/include
	${MPIF90} -Mpreprocess -fPIC -module ${BUILDDIR}/include ${DEFINES} ${INCLUDES} -c -o $@ $<

${CUDECOMPLIB}: ${OBJ}
	@mkdir -p ${BUILDDIR}/lib
	${NVCC} -shared ${NVFLAGS} ${INCLUDES} ${LIBS} -o $@ $^ ${STATIC_LIBS}

${CUDECOMPFLIB}: ${CUDECOMPMOD}
	@mkdir -p ${BUILDDIR}/lib
	${MPIF90} -shared ${INCLUDES} -o $@ $^

clean:
	rm -f ${CUDECOMPLIB} ${CUDECOMPFLIB} ${BUILDDIR}/*.o ${BUILDDIR}/include/*.mod ${BUILDDIR}/include/*.h
	cd ${TESTDIR}; make clean
	cd ${BENCHMARKDIR}; make clean
	cd ${EXAMPLEDIR}; make clean
