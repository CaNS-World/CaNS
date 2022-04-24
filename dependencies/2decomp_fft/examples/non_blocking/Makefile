include ../../src/Makefile.inc

INCLUDE = -I../../include
LIBS = -L../../lib -l2decomp_fft $(LIBFFT)

all: blocking non_blocking

blocking: blocking.o
	$(F90) -o $@ $@.o $(LIBS)

non_blocking: non_blocking.o
	$(F90) -o $@ $@.o $(LIBS)

clean:
	rm -f *.o blocking non_blocking

%.o : %.f90
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) -c $<
