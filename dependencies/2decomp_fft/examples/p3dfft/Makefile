include ../../src/Makefile.inc

P3DFFT_HOME=$(HOME)/software/build/p3dfft-2.7.9-dimsc
FFTW3_HOME=$(HOME)/software/build/fftw-3.3

INCLUDE = -I../../include -I$(P3DFFT_HOME)/include
LIBS = -L../../lib -l2decomp_fft -L$(P3DFFT_HOME)/lib -lp3dfft $(LIBFFT) -L$(FFTW3_HOME)/lib -lfftw3

OBJ = p3dfft.o

p3dfft: $(OBJ)
	$(F90) -o $@ $(OBJ) $(LIBS)

clean:
	rm -f *.o p3dfft

%.o : %.f90
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) -c $<
