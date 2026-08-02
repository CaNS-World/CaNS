#!/usr/bin/env python
#
# SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
# SPDX-License-Identifier: MIT
#
def read_single_field_hdf5(data_dir,filenamei,varname=""):
    import os
    import h5py
    import numpy as np
    #
    # setting up some parameters
    #
    iprecision = 8            # precision of the real-valued data
    r0 = np.array([0.,0.,0.]) # domain origin
    precision  = 'float64'
    if(iprecision == 4): precision = 'float32'
    #
    # read field file
    #
    filepath = os.path.join(data_dir,filenamei)
    hf = h5py.File(filepath,'r')
    if(len(varname) == 0):
        fields = list(hf["fields"].keys())
        if(len(fields) != 1):
            raise ValueError("varname must be provided when the HDF5 file stores multiple fields")
        varname = fields[0]
    data = np.transpose(np.asarray(hf["fields/"+varname]))
    if("meta/lo" in hf):
        lo = np.asarray(hf["meta/lo"],dtype=int)
        hi = np.asarray(hf["meta/hi"],dtype=int)
        nskip = np.asarray(hf["meta/nskip"],dtype=int)
    else:
        ng_local = np.array(data.shape,dtype=int)
        lo = np.array([1,1,1],dtype=int)
        hi = ng_local.copy()
        nskip = np.array([1,1,1],dtype=int)
    hf.close()
    #
    # read geometry file
    #
    geofile  = data_dir+"/geometry.out"
    geo = np.loadtxt(geofile, comments = "!", max_rows = 2)
    ng = geo[0,:].astype('int')
    l  = geo[1,:]
    dl = l/(1.*ng)
    #
    # read and generate grid
    #
    xp = np.linspace(r0[0]+dl[0]/2.,r0[0]+l[0]-dl[0]/2.,ng[0]) # centered grid
    yp = np.linspace(r0[1]+dl[1]/2.,r0[1]+l[1]-dl[1]/2.,ng[1]) # centered grid
    zp = np.linspace(r0[2]+dl[2]/2.,r0[2]+l[2]-dl[2]/2.,ng[2]) # centered grid
    xu = xp + dl[0]/2. # staggered grid
    yv = yp + dl[1]/2. # staggered grid
    zw = zp + dl[2]/2. # staggered grid
    if(os.path.exists(data_dir+"/grid.h5")):
        hf = h5py.File(data_dir+"/grid.h5","r")
        zp = np.asarray(hf["rc"])
        zw = np.asarray(hf["rf"])
        hf.close()
    elif(os.path.exists(data_dir+"/grid.bin")):
        f = open(data_dir+'/grid.bin','rb')
        grid_z = np.fromfile(f,dtype=precision)
        f.close()
        grid_z = np.reshape(grid_z,(ng[2],4),order='F')
        zp = r0[2] + grid_z[:,2] # centered  z grid
        zw = r0[2] + grid_z[:,3] # staggered z grid
    #
    # reshape grid
    #
    xp = xp[lo[0]-1:hi[0]:nskip[0]]
    yp = yp[lo[1]-1:hi[1]:nskip[1]]
    zp = zp[lo[2]-1:hi[2]:nskip[2]]
    xu = xu[lo[0]-1:hi[0]:nskip[0]]
    yv = yv[lo[1]-1:hi[1]:nskip[1]]
    zw = zw[lo[2]-1:hi[2]:nskip[2]]
    return data, xp, yp, zp, xu, yv, zw
if __name__ == "__main__":
    filenamei = input("Name of the HDF5 file written by CaNS (e.g. vex_fld_0000000.h5)]: ")
    varname = input("Name of the HDF5 field dataset []: ") or ""
    data, xp, yp, zp, xu, yv, zw = read_single_field_hdf5("./",filenamei,varname)
