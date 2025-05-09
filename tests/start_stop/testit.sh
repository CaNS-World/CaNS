#!/bin/bash
#
set -e # exit immediately if any command fails
#
TESTDIR=$(pwd)
CANSDIR=$(pwd)/../..
SRCDIR=$CANSDIR/src
RUNDIR=$CANSDIR/run
UTILSDIR=$CANSDIR/utils

MPIRUN_OPTIONS='--oversubscribe'
if mpirun --version 2>&1 | grep -qi intel; then MPIRUN_OPTIONS=""; fi
MPIRUN="mpirun -n 4 $MPIRUN_OPTIONS"

cd $RUNDIR
mkdir -p data
#
cp $TESTDIR/input-*.nml $RUNDIR

NAME_ONE="one-step"
NAME_TWO_1="two-step-1"
NAME_TWO_2="two-step-2"

for name in $NAME_ONE $NAME_TWO_1 $NAME_TWO_2; do
  echo "INFO: Running CaNS, case $name"
  mv input-${name}.nml input.nml
  ${MPIRUN} ./cans 1> log_file.log 2> err_log.log
  (head -n 50 log_file.log; echo -e "\n[...output omitted...]\n"; tail -n 50 log_file.log) # report first and last N lines of the log file
  dirname="data-$(echo $name | sed 's/-[1-2]//g')/"
  mkdir -p $dirname
  if [[ $name == *"two-step"* ]]; then
    cp -r data/* $dirname
  else
    mv data/* $dirname
  fi
done

cp $TESTDIR/*.* . && cp $UTILSDIR/read_binary_data/python/read_single_field_binary.py .
echo "INFO: Running comparison to reference"
pytest test.py
rm -rf $RUNDIR/data/*.* data-*
