#!/bin/bash
#
# 1) generate grid files
#
gfortran gen_grid.f90 -o a.out && ./a.out && rm -rf a.out
#
# 2) generate xdmf file
#
gfortran gen_xdmf.f90 -o a.out && ./a.out && rm -rf a.out
