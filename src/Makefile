# src directory
SRCDIR := $(pwd)
RUNDIR := $(SRCDIR)'../run'
#
# fortran compiler
FC   = mpif90
# flags
BIG := # -mcmodel=medium
DBG :=  -g -fbacktrace -Wall -O0 -Wextra -pedantic -fcheck=all
OMP := # -fopenmp
OPT := #-O3
OTH := -DDEBUG -DTIMING #-DIMPDIFF
PROF := #-pg
FFLAGS := $(OPT) $(BIG) $(DBG) $(PROF) $(OMP) $(OTH)
TARGET = cans
#
OBJ = bound.o chkdiv.o chkdt.o common_mpi.o correc.o debug.o fft.o fftw.o fillps.o initflow.o initgrid.o initmpi.o initsolver.o load.o main.o mom.o momd.o moms.o output.o param.o rk.o sanity.o solver.o decomp_2d.o io.o
#
LIBS = -L /usr/local/Cellar/fftw/3.3.8_1/lib  -lfftw3 -lfftw3_threads
#
all: $(TARGET)
#
$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) $(LIBS) -o $(TARGET)
#
.PHONY: clean
clean:
	rm -rf *.o *.mod *dSYM $(TARGET)
# dependencies
bound.o: bound.f90 common_mpi.o
	$(FC) $(FFLAGS) -cpp -c  $<
chkdiv.o: chkdiv.f90 common_mpi.o
	$(FC) $(FFLAGS) -cpp -c  $<
chkdt.o: chkdt.f90 common_mpi.o
	$(FC) $(FFLAGS) -cpp -c  $<
common_mpi.o: common_mpi.f90 param.o
	$(FC) $(FFLAGS) -cpp -c  $<
correc.o: correc.f90
	$(FC) $(FFLAGS) -cpp -c  $<
debug.o: debug.f90 common_mpi.o
	$(FC) $(FFLAGS) -cpp -c  $<
fillps.o: fillps.f90
	$(FC) $(FFLAGS) -cpp -c  $<
fft.o: fft.f90 common_mpi.o fftw.o
	$(FC) $(FFLAGS) -cpp -c  $<
fftw.o: fftw.f90
	$(FC) $(FFLAGS) -cpp -c  $<
initflow.o: initflow.f90 common_mpi.o decomp_2d.o
	$(FC) $(FFLAGS) -cpp -c  $<
initgrid.o: initgrid.f90 param.o
	$(FC) $(FFLAGS) -cpp -c  $<
initmpi.o: initmpi.f90 decomp_2d.o common_mpi.o
	$(FC) $(FFLAGS) -cpp -c  $<
initsolver.o: initsolver.f90 common_mpi.o fft.o param.o
	$(FC) $(FFLAGS) -cpp -c  $<
load.o: load.f90 common_mpi.o decomp_2d.o io.o
	$(FC) $(FFLAGS) -cpp -c  $<
main.o: main.f90 bound.o chkdiv.o chkdt.o debug.o fft.o fillps.o initflow.o initflow.o initgrid.o initmpi.o initsolver.o load.o output.o param.o rk.o sanity.o solver.o out1d.h90 out2d.h90 out3d.h90
	$(FC) $(FFLAGS) -cpp -c  $<
mom.o: mom.f90 common_mpi.o
	$(FC) $(FFLAGS) -cpp -c  $<
momd.o: momd.f90 common_mpi.o
	$(FC) $(FFLAGS) -cpp -c  $<
moms.o: moms.f90
	$(FC) $(FFLAGS) -cpp -c  $<
output.o: output.f90 io.o
	$(FC) $(FFLAGS) -cpp -c  $<
param.o: param.f90 setup.h90 bc.h90
	$(FC) $(FFLAGS) -cpp -c  $<
rk.o: rk.f90 debug.o mom.o momd.o moms.o param.o
	$(FC) $(FFLAGS) -cpp -c  $<
sanity.o: sanity.f90 bound.o chkdiv.o common_mpi.o correc.o debug.o decomp_2d.o fft.o fillps.o initflow.o initmpi.o initsolver.o param.o solver.o
	$(FC) $(FFLAGS) -cpp -c  $<
solver.o: solver.f90 fft.o param.o
	$(FC) $(FFLAGS) -cpp -c  $<
# 2decomp&fft library
FFLAGS_2DECOMP := -DDOUBLE_PREC# -DOVERWRITE -DEVEN
2decompdir=2decomp/
decomp_2d.o: $(2decompdir)decomp_2d.f90
	$(FC) $(FFLAGS) -cpp $(FFLAGS_2DECOMP) -c $<
io.o: $(2decompdir)io.f90 decomp_2d.o
	$(FC) $(FFLAGS) -cpp $(FFLAGS_2DECOMP) -c $<
