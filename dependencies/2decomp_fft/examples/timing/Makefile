include ../../src/Makefile.inc

INCLUDE = -I../../include
LIBS = -L../../lib -l2decomp_fft $(LIBFFT)

OBJ = timing.o

timing: $(OBJ)
	$(F90) -o $@ $(OBJ) $(LIBS)

clean:
	rm -f *.o timing

%.o : %.f90
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) -c $<
