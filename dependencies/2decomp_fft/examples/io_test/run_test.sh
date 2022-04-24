#!/bin/sh

make
if [ $? -ne 0 ] ; then
    echo "================================================="
    echo "Failed to build the applications. Fix them first!"
    echo "================================================="
    exit 1
fi

echo " "
echo "writing data files using MPI-IO..."
echo "mpirun -np 12 ./io_test"
mpirun -np 12 ./io_test

echo " "
echo "reading data files back (different number of processes)..."
echo "mpirun -np 6 ./io_read"
mpirun -np 6 ./io_read

# The files written by MPI-IO should be independent on # processes
echo " "
echo "*** testing write_var..."
echo "mpirun -np 20 ./io_var_test 5 4"
mpirun -np 20 ./io_var_test 5 4
echo "mpirun -np 12 ./io_var_test 4 3"
mpirun -np 12 ./io_var_test 4 3
echo "mpirun -np 6 ./io_var_test 3 2"
mpirun -np 6 ./io_var_test 3 2
echo "mpirun -np 2 ./io_var_test 2 1"
mpirun -np 2 ./io_var_test 2 1
diff -s io_var_data.020 io_var_data.012 
diff -s io_var_data.020 io_var_data.006
diff -s io_var_data.020 io_var_data.002

echo " "
echo "*** testing write_plane..."
echo "mpirun -np 12 ./io_plane_test"
mpirun -np 12 ./io_plane_test
diff -s x_pencil-x_plane.dat y_pencil-x_plane.dat
diff -s x_pencil-x_plane.dat z_pencil-x_plane.dat
diff -s x_pencil-y_plane.dat y_pencil-y_plane.dat
diff -s x_pencil-y_plane.dat z_pencil-y_plane.dat
diff -s x_pencil-z_plane.dat y_pencil-z_plane.dat
diff -s x_pencil-z_plane.dat z_pencil-z_plane.dat

echo " "
echo "Tests PASSED, unless errors reported"

