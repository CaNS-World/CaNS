#!/usr/bin/env python
#
# SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
# SPDX-License-Identifier: MIT
#
def read_restart_file_hdf5(filenamei):
    import h5py
    import numpy as np
    #
    # read checkpoint HDF5 file
    #
    with h5py.File(filenamei, "r") as hf:
        u = np.asarray(hf["fields/u"])
        v = np.asarray(hf["fields/v"])
        w = np.asarray(hf["fields/w"])
        p = np.asarray(hf["fields/p"])
        time = float(np.asarray(hf["meta/time"])[0])
        istep = int(np.asarray(hf["meta/istep"])[0])
    return u, v, w, p, time, istep
if __name__ == "__main__":
    filenamei = input("Name of the HDF5 restart file written by CaNS [fld.h5]: ") or "fld.h5"
    u, v, w, p, time, istep = read_restart_file_hdf5(filenamei)