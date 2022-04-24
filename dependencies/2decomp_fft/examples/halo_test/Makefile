include ../../src/Makefile.inc

INCLUDE = -I../../include
LIBS = -L../../lib -l2decomp_fft

OBJ = halo_test.o

halo_test: $(OBJ)
	$(F90) -o $@ $(OBJ) $(LIBS)

clean:
	rm -f *.o halo_test

%.o : %.f90
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) -c $<
