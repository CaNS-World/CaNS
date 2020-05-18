#!/usr/bin/env python
import numpy as np
#
# setting up some parameters
#
iprecision = 8            # precision of the real-valued data
r0 = np.array([0.,0.,0.]) # domain origin
non_uniform_grid = True
#
# read geometry file
#
geofile  = "geometry.out"
geo = np.loadtxt(geofile, comments = "!", max_rows = 2)
ng = geo[0,:].astype('int')
l  = geo[1,:]
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
    if(    iprecision == 4):
        grid_z = np.fromfile(f,dtype='float32')
    else:
        grid_z = np.fromfile(f,dtype='float64')
    f.close()
    grid_z = np.reshape(grid_z,(ng[2],4),order='F')
    zp = r0[2] + grid_z[:,2] # centered  z grid
    zw = r0[2] + grid_z[:,1] # staggered z grid
#
# read binary file
#
filenamei   = input("Name of the binary file written by CaNS (e.g. vex_fld_0000000.bin)]: ")
iskipx      = input("Data saved every (ix, iy, iz) points. Value of ix? [1]: ") or "1"
iskipy      = input("Data saved every (ix, iy, iz) points. Value of iy? [1]: ") or "1"
iskipz      = input("Data saved every (ix, iy, iz) points. Value of iz? [1]: ") or "1"
iskip       = np.array([iskipx,iskipy,iskipz]).astype(int)
n           = (ng[:]/iskip[:]).astype(int)
data        = np.zeros([n[0],n[1],n[2]]) # u[:,:,:],v[:,:,:],w[:,:,:],p[:,:,:]
fld         = np.fromfile(filenamei,dtype='float64')
data[:,:,:] = np.reshape(fld,(n[0],n[1],n[2]),order='F')
#
# reshape grid
#
xp = xp[0:ng[0]:iskip[0]]
yp = yp[0:ng[1]:iskip[1]]
zp = zp[0:ng[2]:iskip[2]]
xu = xu[0:ng[0]:iskip[0]]
yv = yv[0:ng[1]:iskip[1]]
zw = zw[0:ng[2]:iskip[2]]
