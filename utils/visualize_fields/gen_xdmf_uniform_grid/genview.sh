#!/bin/bash
#
# generate xdmf file
#
gfortran gen_xdmf.f90 -o a.out && ./a.out && rm -rf a.out
