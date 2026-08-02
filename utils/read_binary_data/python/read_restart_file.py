# -
#
# SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
# SPDX-License-Identifier: MIT
#
# -
#!/usr/bin/env python
def read_restart_file(filenamei):
    import os
    import numpy as np
    #
    # setting up some parameters
    #
    iprecision = 8            # precision of the real-valued data
    r0 = np.array([0.,0.,0.]) # domain origin
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
    xp = np.linspace(r0[0]+dl[0]/2.,r0[0]+l[0]-dl[0]/2.,ng[0]) # centered grid
    yp = np.linspace(r0[1]+dl[1]/2.,r0[1]+l[1]-dl[1]/2.,ng[1]) # centered grid
    zp = np.linspace(r0[2]+dl[2]/2.,r0[2]+l[2]-dl[2]/2.,ng[2]) # centered grid
    xu = xp + dl[0]/2. # staggered grid
    yv = yp + dl[1]/2. # staggered grid
    zw = zp + dl[2]/2. # staggered grid
    if(os.path.exists('grid.bin')):
        f   = open('grid.bin','rb')
        grid_z = np.fromfile(f,dtype=precision)
        f.close()
        grid_z = np.reshape(grid_z,(ng[2],4),order='F')
        zp = r0[2] + grid_z[:,2] # centered  z grid
        zw = r0[2] + grid_z[:,3] # staggered z grid
    #
    # read checkpoint binary file
    #
    disp = np.prod(ng)
    fldinfo = np.zeros([2])
    legacy_size = iprecision*(4*disp+2)
    def split_prefix(filename):
        root, ext = os.path.splitext(filename)
        if(ext == ".bin"):
            for suffix in ["_u","_v","_w","_p"]:
                if(root.endswith(suffix)):
                    return root[:-len(suffix)]
            return root
        return filename
    def read_one(filename):
        with open(filename,'rb') as f:
            fld = np.fromfile(f,dtype=precision,count=disp)
            info = np.fromfile(f,dtype=precision,count=2)
        return np.reshape(fld,(ng[0],ng[1],ng[2]),order='F'), info
    if(os.path.exists(filenamei) and os.path.getsize(filenamei) == legacy_size):
        #
        # u, v, w, p
        #
        data = np.zeros([ng[0],ng[1],ng[2],4])
        offset = 0
        with open(filenamei,'rb') as f:
            for q in range(4):
                f.seek(offset)
                fld = np.fromfile(f,dtype=precision,count=disp)
                data[:,:,:,q] = np.reshape(fld,(ng[0],ng[1],ng[2]),order='F')
                offset += iprecision*disp
            f.seek(offset)
            fldinfo[:] = np.fromfile(f,dtype=precision,count=2)
        u = data[:,:,:,0]
        v = data[:,:,:,1]
        w = data[:,:,:,2]
        p = data[:,:,:,3]
    else:
        prefix = split_prefix(filenamei)
        u, fldinfo = read_one(prefix+"_u.bin")
        v, _       = read_one(prefix+"_v.bin")
        w, _       = read_one(prefix+"_w.bin")
        p, _       = read_one(prefix+"_p.bin")
    time  =     fldinfo[0]
    istep = int(fldinfo[1])
    return u, v, w, p, time, istep

if __name__ == "__main__":
    filenamei = input("Name of the binary restart file written by CaNS [fld.bin]: ") or "fld.bin"
    u, v, w, p, time, istep = read_restart_file(filenamei)
