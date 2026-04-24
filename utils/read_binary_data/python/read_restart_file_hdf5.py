#!/usr/bin/env python
#
# SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
# SPDX-License-Identifier: MIT
#
def read_restart_file_hdf5(filenamei):
    import os
    import h5py
    import numpy as np
    #
    # read checkpoint HDF5 file
    #
    def split_prefix(filename):
        root, ext = os.path.splitext(filename)
        if(ext == ".h5"):
            for suffix in ["_u","_v","_w","_p"]:
                if(root.endswith(suffix)):
                    return root[:-len(suffix)]
            return root
        return filename
    def read_one(filename,fieldname):
        with h5py.File(filename, "r") as hf:
            fld = np.transpose(np.asarray(hf["fields/"+fieldname]))
            time = float(np.asarray(hf["meta/time"])[0])
            istep = int(np.asarray(hf["meta/istep"])[0])
        return fld, time, istep
    with h5py.File(filenamei, "r") as hf:
        fieldnames = list(hf["fields"].keys())
    if(all(fieldname in fieldnames for fieldname in ["u","v","w","p"])):
        with h5py.File(filenamei, "r") as hf:
            u = np.transpose(np.asarray(hf["fields/u"]))
            v = np.transpose(np.asarray(hf["fields/v"]))
            w = np.transpose(np.asarray(hf["fields/w"]))
            p = np.transpose(np.asarray(hf["fields/p"]))
            time = float(np.asarray(hf["meta/time"])[0])
            istep = int(np.asarray(hf["meta/istep"])[0])
    else:
        prefix = split_prefix(filenamei)
        u, time, istep = read_one(prefix+"_u.h5","u")
        v, _   , _     = read_one(prefix+"_v.h5","v")
        w, _   , _     = read_one(prefix+"_w.h5","w")
        p, _   , _     = read_one(prefix+"_p.h5","p")
    return u, v, w, p, time, istep
if __name__ == "__main__":
    filenamei = input("Name of the HDF5 restart file written by CaNS [fld.h5]: ") or "fld.h5"
    u, v, w, p, time, istep = read_restart_file_hdf5(filenamei)
