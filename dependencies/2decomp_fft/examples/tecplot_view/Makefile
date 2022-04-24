include ../../src/Makefile.inc

INCLUDE = -I../../include
LIBS = -L../../lib -l2decomp_fft

OBJ = tecplot_view.o

tecplot_view: $(OBJ)
	$(F90) -o $@ $(OBJ) $(LIBS)

clean:
	rm -f *.o tecplot_view

%.o : %.f90
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) -c $<
