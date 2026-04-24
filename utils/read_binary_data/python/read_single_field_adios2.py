#!/usr/bin/env python
#
# SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
# SPDX-License-Identifier: MIT
#
def read_single_field_adios2(data_dir,filenamei,varname=""):
    import os
    os.environ.setdefault("ADIOS2_ALWAYS_USE_MPI","1")
    import adios2
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
    FileReader = getattr(adios2,"FileReader",None)
    if(FileReader):
        with FileReader(filepath) as fh:
            available_variables = fh.available_variables()
            if(len(varname) == 0):
                meta_names = {"time","istep","lo","hi","nskip","x","y","z"}
                fields = [name for name in available_variables.keys() if name not in meta_names]
                if(len(fields) != 1):
                    raise ValueError("varname must be provided when the ADIOS2 file stores multiple fields")
                varname = fields[0]
            data = np.transpose(np.asarray(fh.read(varname)))
            if("lo" in available_variables):
                lo = np.asarray(fh.read("lo"),dtype=int)
                hi = np.asarray(fh.read("hi"),dtype=int)
                nskip = np.asarray(fh.read("nskip"),dtype=int)
            else:
                ng_local = np.array(data.shape,dtype=int)
                lo = np.array([1,1,1],dtype=int)
                hi = ng_local.copy()
                nskip = np.array([1,1,1],dtype=int)
    else:
        adios = adios2.ADIOS()
        io = adios.DeclareIO("reader")
        io.SetEngine("BP5")
        engine = io.Open(filepath,adios2.Mode.Read)
        engine.BeginStep()
        available_variables = io.AvailableVariables()
        if(len(varname) == 0):
            meta_names = {"time","istep","lo","hi","nskip","x","y","z"}
            fields = [name for name in available_variables.keys() if name not in meta_names]
            if(len(fields) != 1):
                raise ValueError("varname must be provided when the ADIOS2 file stores multiple fields")
            varname = fields[0]
        def read_bp_var(name,dtype=None):
            var = io.InquireVariable(name)
            if(var is None):
                raise KeyError(name)
            shape = tuple(var.Shape())
            if(dtype is None):
                dtype = np.float32 if var.Type() == "float" else np.float64
                if("int" in var.Type()): dtype = np.int32
            arr = np.empty(shape,dtype=dtype)
            engine.Get(var,arr,adios2.Mode.Sync)
            return arr
        data = np.transpose(read_bp_var(varname))
        if("lo" in available_variables):
            lo = read_bp_var("lo",np.int32).astype(int)
            hi = read_bp_var("hi",np.int32).astype(int)
            nskip = read_bp_var("nskip",np.int32).astype(int)
        else:
            ng_local = np.array(data.shape,dtype=int)
            lo = np.array([1,1,1],dtype=int)
            hi = ng_local.copy()
            nskip = np.array([1,1,1],dtype=int)
        engine.EndStep()
        engine.Close()
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
        import h5py
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
    filenamei = input("Name of the ADIOS2 file written by CaNS (e.g. vex_fld_0000000.bp)]: ")
    varname = input("Name of the ADIOS2 field variable []: ") or ""
    data, xp, yp, zp, xu, yv, zw = read_single_field_adios2("./",filenamei,varname)
