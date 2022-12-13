# -
#
# SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
# SPDX-License-Identifier: MIT
#
# -
#!/usr/bin/env python
def read_restart_file(filenamei):
    import numpy as np
    #
    # setting up some parameters
    #
    iprecision = 8            # precision of the real-valued data
    r0 = np.array([0.,0.,0.]) # domain origin
    non_uniform_grid = True
    precision  = 'float64'
    if(iprecision == 4): precision = 'float32'
    #
    # read geometry file
    #
    geofile  = "geometry.out"
    data = np.loadtxt(geofile, comments = "!", max_rows = 2)
    ng = data[0,:].astype('int')
    l  = data[1,:]
    dl = l/(1.*ng)
    #
    # read and generate grid
    #
    xp = np.arange(r0[0]+dl[0]/2.,r0[0]+l[0],dl[0]) # centered  x grid
    yp = np.arange(r0[1]+dl[1]/2.,r0[1]+l[1],dl[1]) # centered  y grid
    zp = np.arange(r0[2]+dl[2]/2.,r0[2]+l[2],dl[2]) # centered  z grid
    xu = xp + dl[0]/2.                              # staggered x grid
    yv = yp + dl[1]/2.                              # staggered y grid
    zw = zp + dl[2]/2.                              # staggered z grid
    if(non_uniform_grid):
        f   = open('grid.bin','rb')
        grid_z = np.fromfile(f,dtype=precision)
        f.close()
        grid_z = np.reshape(grid_z,(ng[2],4),order='F')
        zp = r0[2] + grid_z[:,2] # centered  z grid
        zw = r0[2] + grid_z[:,3] # staggered z grid
    #
    # read checkpoint binary file
    #
    offset     = 0
    disp       = np.prod(ng)
    data       = np.zeros([ng[0],ng[1],ng[2],4]) # u[:,:,:],v[:,:,:],w[:,:,:],p[:,:,:]
    fldinfo    = np.zeros([2])
    with open(filenamei,'rb') as f:
        for p in range(4):
            f.seek(offset)
            fld = np.fromfile(f,dtype=precision,count=disp)
            data[:,:,:,p] = np.reshape(fld,(ng[0],ng[1],ng[2]),order='F')
            offset += iprecision*disp
        f.seek(offset)
        fldinfo[:] = np.fromfile(f,dtype=precision,count=2)
    f.close
    #
    # store data in arrays
    #
    u = data[:,:,:,0]
    v = data[:,:,:,1]
    w = data[:,:,:,2]
    p = data[:,:,:,3]
    time  =     fldinfo[0]
    istep = int(fldinfo[1])
    return u, v, w, p, time, istep

if __name__ == "__main__":
    filenamei = input("Name of the binary restart file written by CaNS [fld.bin]: ") or "fld.bin"
    u, v, w, p, time, istep = read_restart_file(filenamei)
