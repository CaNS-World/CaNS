#!/bin/bash
TESTDIR=$(pwd)
CANSDIR=$(pwd)/../..
SRCDIR=$CANSDIR/src
RUNDIR=$CANSDIR/run
UTILSDIR=$CANSDIR/utils
rm -rf $RUNDIR
echo "Compiling ..."
sleep 2
cp $TESTDIR/Makefile $SRCDIR && cd $SRCDIR && make clean && make -j run
cp $TESTDIR/dns.in $RUNDIR && cd $RUNDIR
echo "Running CaNS..."
sleep 2
mpirun -n 4 ./cans
cp $TESTDIR/*.* data/ && cp $UTILSDIR/read_binary_data/python/read_single_field_binary.py $RUNDIR/data/ && cd $RUNDIR/data/
echo "Running test..."
sleep 2
pytest test.py
rm -rf ../../run
