#!/usr/bin/env python
#
# SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
# SPDX-License-Identifier: MIT
#
def read_restart_file_adios2(filenamei):
    import adios2
    import numpy as np
    #
    # read checkpoint ADIOS2 file
    #
    with adios2.FileReader(filenamei) as fh:
        u = np.asarray(fh.read("u"))
        v = np.asarray(fh.read("v"))
        w = np.asarray(fh.read("w"))
        p = np.asarray(fh.read("p"))
        time = float(np.asarray(fh.read("time"))[0])
        istep = int(np.asarray(fh.read("istep"))[0])
    return u, v, w, p, time, istep
if __name__ == "__main__":
    filenamei = input("Name of the ADIOS2 restart file written by CaNS [fld.bp]: ") or "fld.bp"
    u, v, w, p, time, istep = read_restart_file_adios2(filenamei)
