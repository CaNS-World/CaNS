# -
#
# SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
# SPDX-License-Identifier: MIT
#
# -
def estimate_memory_footprint(nx,ny,nz,nproc=1,rp=8,gp=8,is_impdiff=False,is_z_periodic=False,is_transpose_in_place=False):
    n = nx*ny*nz/nproc
    footprint = 0
    footprint += n*5*rp # main arrays
    footprint += n*6*rp # prediction velocity arrays
    footprint += n*2*gp # solver pressure buffer arrays
    footprint += n*2*gp # other work arrays
    if(is_impdiff):
        footprint += n*3*rp # implicit diffusion extra arrays
    if(is_z_periodic):
        footprint += n*2*gp # extra arrays for z peridicity
    if(is_transpose_in_place):
        footprint -= n*1*gp # less memory footprint
    return footprint

if __name__ == '__main__':
    #
    nx         = int(input("Number of grid points along x [default: 256]? ") or 256)
    ny         = int(input("Number of grid points along y [default: 256]? ") or 256)
    nz         = int(input("Number of grid points along z [default: 256]? ") or 256)
    nproc      = int(input("Number of GPUs [default: 1]? ") or 1)
    gp         = int(input("Poisson solver precision (8: double or 4: single) [default: 8]? ") or 8)
    is_impdiff = int(input("Implicit Diffusion (0: no, else: yes) [default: 0]? ") or 0) != 0
    #
    footprint = estimate_memory_footprint(nx,ny,nz,nproc,rp=8,gp=gp,is_impdiff=is_impdiff)
    print("Estimated memory footprint (Gb): ", footprint/1024.**3)
    print("Estimated memory footprint (bytes per grid point): ", footprint/(nx*ny*nz/nproc))
