# fortran compiler
FC   = mpif90
#
BIG := # -mcmodel=medium
DBG := # -g -fbacktrace -Wall -O0 -Wextra -pedantic -fcheck=all
OMP := # -fopenmp
OPT := -O3
OTH := -DDEBUG -DTIMING #-DIMPDIFF
PROF := #-pg
FFLAGS := -cpp $(OPT) $(BIG) $(DBG) $(PROF) $(OMP) $(OTH)

#
TARGET = cans
#
OBJ = bound.o chkdiv.o chkdt.o common_mpi.o correc.o debug.o fft.o fftw.o fillps.o initflow.o initgrid.o initmpi.o initsolver.o load.o main.o mom.o momd.o moms.o output.o param.o rk.o sanity.o solver.o decomp_2d.o io.o
#
LIBS = -lfftw3 -lfftw3_threads -llapack
#
all: $(TARGET)
#
$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) $(LIBS) -o $(TARGET)
#
.PHONY: clean
clean:
	rm -rf *.o *.mod *dSYM $(TARGET)
#
bound.o: bound.f90 common_mpi.o
	$(FC) $(FFLAGS) -c  bound.f90 
chkdiv.o: chkdiv.f90 common_mpi.o
	$(FC) $(FFLAGS) -c  chkdiv.f90 
chkdt.o: chkdt.f90 common_mpi.o
	$(FC) $(FFLAGS) -c  chkdt.f90 
common_mpi.o: common_mpi.f90 param.o
	$(FC) $(FFLAGS) -c  common_mpi.f90 
correc.o: correc.f90
	$(FC) $(FFLAGS) -c  correc.f90 
debug.o: debug.f90 common_mpi.o
	$(FC) $(FFLAGS) -c  debug.f90 
fillps.o: fillps.f90
	$(FC) $(FFLAGS) -c  fillps.f90
fft.o: fft.f90 common_mpi.o fftw.o
	$(FC) $(FFLAGS) -c  fft.f90
fftw.o: fftw.f90
	$(FC) $(FFLAGS) -c  fftw.f90
initflow.o: initflow.f90 common_mpi.o decomp_2d.o
	$(FC) $(FFLAGS) -c  initflow.f90 
initgrid.o: initgrid.f90 param.o
	$(FC) $(FFLAGS) -c  initgrid.f90 
initmpi.o: initmpi.f90 decomp_2d.o common_mpi.o
	$(FC) $(FFLAGS) -c  initmpi.f90
initsolver.o: initsolver.f90 common_mpi.o fft.o param.o
	$(FC) $(FFLAGS) -c  initsolver.f90 
load.o: load.f90 common_mpi.o decomp_2d.o io.o
	$(FC) $(FFLAGS) -c  load.f90 
main.o: main.f90 bound.o chkdiv.o chkdt.o debug.o fft.o fillps.o initflow.o initflow.o initgrid.o initmpi.o initsolver.o load.o output.o param.o rk.o sanity.o solver.o
	$(FC) $(FFLAGS) -c  main.f90
mom.o: mom.f90 common_mpi.o
	$(FC) $(FFLAGS) -c  mom.f90
momd.o: momd.f90 common_mpi.o
	$(FC) $(FFLAGS) -c  momd.f90
moms.o: moms.f90
	$(FC) $(FFLAGS) -c  moms.f90
output.o: output.f90 io.o
	$(FC) $(FFLAGS) -c  output.f90
param.o: param.f90 setup.h90 bc.h90
	$(FC) $(FFLAGS) -c  param.f90
rk.o: rk.f90 debug.o mom.o momd.o moms.o param.o
	$(FC) $(FFLAGS) -c  rk.f90
sanity.o: sanity.f90 bound.o chkdiv.o common_mpi.o correc.o debug.o decomp_2d.o fft.o fillps.o initflow.o initmpi.o initsolver.o param.o solver.o
	$(FC) $(FFLAGS) -c  sanity.f90
solver.o: solver.f90 fft.o param.o
	$(FC) $(FFLAGS) -c  solver.f90
#
FFLAGS_2DECOMP := -DDOUBLE_PREC# -DOVERWRITE -DEVEN
2decompdir=2decomp/
decomp_2d.o: $(2decompdir)decomp_2d.f90
	$(FC) $(FFLAGS) $(FFLAGS_2DECOMP) -c $(2decompdir)decomp_2d.f90
io.o: $(2decompdir)io.f90 decomp_2d.o
	$(FC) $(FFLAGS) $(FFLAGS_2DECOMP) -c $(2decompdir)io.f90
