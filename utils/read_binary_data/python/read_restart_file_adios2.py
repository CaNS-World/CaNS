#!/usr/bin/env python
#
# SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
# SPDX-License-Identifier: MIT
#
def read_restart_file_adios2(filenamei):
    import os
    os.environ.setdefault("ADIOS2_ALWAYS_USE_MPI","1")
    import adios2
    import numpy as np
    #
    # read checkpoint ADIOS2 file
    #
    def split_prefix(filename):
        root, ext = os.path.splitext(filename)
        if(ext == ".bp"):
            for suffix in ["_u","_v","_w","_p"]:
                if(root.endswith(suffix)):
                    return root[:-len(suffix)]
            return root
        return filename
    FileReader = getattr(adios2,"FileReader",None)
    if(FileReader):
        def available(filename):
            with FileReader(filename) as fh:
                return fh.available_variables()
        def read_one(filename,fieldname):
            with FileReader(filename) as fh:
                fld = np.transpose(np.asarray(fh.read(fieldname)))
                time = float(np.asarray(fh.read("time"))[0])
                istep = int(np.asarray(fh.read("istep"))[0])
            return fld, time, istep
        available_variables = available(filenamei)
        if(all(fieldname in available_variables for fieldname in ["u","v","w","p"])):
            with FileReader(filenamei) as fh:
                u = np.transpose(np.asarray(fh.read("u")))
                v = np.transpose(np.asarray(fh.read("v")))
                w = np.transpose(np.asarray(fh.read("w")))
                p = np.transpose(np.asarray(fh.read("p")))
                time = float(np.asarray(fh.read("time"))[0])
                istep = int(np.asarray(fh.read("istep"))[0])
        else:
            prefix = split_prefix(filenamei)
            u, time, istep = read_one(prefix+"_u.bp","u")
            v, _   , _     = read_one(prefix+"_v.bp","v")
            w, _   , _     = read_one(prefix+"_w.bp","w")
            p, _   , _     = read_one(prefix+"_p.bp","p")
    else:
        adios = adios2.ADIOS()
        io_count = [0]
        def open_bp(filename):
            io_count[0] += 1
            io = adios.DeclareIO("reader_"+str(io_count[0]))
            io.SetEngine("BP5")
            engine = io.Open(filename,adios2.Mode.Read)
            engine.BeginStep()
            return io,engine
        def read_bp_var(io,engine,name,dtype=None):
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
        def read_one(filename,fieldname):
            io,engine = open_bp(filename)
            fld = np.transpose(read_bp_var(io,engine,fieldname))
            time = float(read_bp_var(io,engine,"time")[0])
            istep = int(read_bp_var(io,engine,"istep",np.int32)[0])
            engine.EndStep()
            engine.Close()
            return fld, time, istep
        io,engine = open_bp(filenamei)
        available_variables = io.AvailableVariables()
        engine.EndStep()
        engine.Close()
        if(all(fieldname in available_variables for fieldname in ["u","v","w","p"])):
            io,engine = open_bp(filenamei)
            u = np.transpose(read_bp_var(io,engine,"u"))
            v = np.transpose(read_bp_var(io,engine,"v"))
            w = np.transpose(read_bp_var(io,engine,"w"))
            p = np.transpose(read_bp_var(io,engine,"p"))
            time = float(read_bp_var(io,engine,"time")[0])
            istep = int(read_bp_var(io,engine,"istep",np.int32)[0])
            engine.EndStep()
            engine.Close()
        else:
            prefix = split_prefix(filenamei)
            u, time, istep = read_one(prefix+"_u.bp","u")
            v, _   , _     = read_one(prefix+"_v.bp","v")
            w, _   , _     = read_one(prefix+"_w.bp","w")
            p, _   , _     = read_one(prefix+"_p.bp","p")
    return u, v, w, p, time, istep
if __name__ == "__main__":
    filenamei = input("Name of the ADIOS2 restart file written by CaNS [fld.bp]: ") or "fld.bp"
    u, v, w, p, time, istep = read_restart_file_adios2(filenamei)
