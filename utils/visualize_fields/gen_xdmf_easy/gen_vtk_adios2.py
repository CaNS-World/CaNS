# -
#
# SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
# SPDX-License-Identifier: MIT
#
# -
#
# this python script can be used to populate CaNS ADIOS2 visualization folders
# with vtk.xml files understood by ParaView/VTK
#
import os
os.environ.setdefault("ADIOS2_ALWAYS_USE_MPI","1")
from pathlib import Path
import numpy as np
import adios2
from xml.etree import ElementTree
from xml.dom import minidom
from xml.etree.ElementTree import Element, SubElement
#
# define some custom parameters, not defined in the DNS code
#
tol_uniform = 1.e-5
meta_names = {"time","istep","lo","hi","nskip","x","y","z"}
#
# helper routines
#
def get_unique_files(logfile):
    dtype_saves = np.dtype([                                                   \
                            ('file' , 'U100'), ('variable', 'U100'),           \
                            ('imin' , int)   , ('jmin' , int), ('kmin' , int), \
                            ('imax' , int)   , ('jmax' , int), ('kmax' , int), \
                            ('istep', int)   , ('jstep', int), ('kstep', int), \
                            ('time', float)  , ('isave', int)                  \
                           ])
    saves = np.loadtxt(logfile,dtype=dtype_saves)
    saves = np.atleast_1d(saves)
    saves = np.unique(saves)
    saves = np.sort(saves,order=['isave','variable'])
    return np.unique(saves['file'])
def get_bp_data(bpdir):
    FileReader = getattr(adios2,"FileReader",None)
    if(FileReader):
        with FileReader(str(bpdir)) as fh:
            available_variables = fh.available_variables()
            fieldnames = [name for name in available_variables.keys() if name not in meta_names]
            if("x" not in available_variables or "y" not in available_variables or "z" not in available_variables):
                raise ValueError("missing x/y/z arrays in "+str(bpdir))
            x = np.asarray(fh.read("x"),dtype=float)
            y = np.asarray(fh.read("y"),dtype=float)
            z = np.asarray(fh.read("z"),dtype=float)
    else:
        adios = adios2.ADIOS()
        io = adios.DeclareIO("reader")
        io.SetEngine("BP5")
        engine = io.Open(str(bpdir),adios2.Mode.Read)
        engine.BeginStep()
        available_variables = io.AvailableVariables()
        fieldnames = [name for name in available_variables.keys() if name not in meta_names]
        if("x" not in available_variables or "y" not in available_variables or "z" not in available_variables):
            raise ValueError("missing x/y/z arrays in "+str(bpdir))
        def read_1d(name):
            var = io.InquireVariable(name)
            arr = np.empty(tuple(var.Shape()),dtype=float)
            engine.Get(var,arr,adios2.Mode.Sync)
            return arr
        x = read_1d("x")
        y = read_1d("y")
        z = read_1d("z")
        engine.EndStep()
        engine.Close()
    return fieldnames,x,y,z
def infer_origin_spacing(arr,name):
    arr = np.asarray(arr,dtype=float)
    if(arr.size == 0):
        raise ValueError("empty coordinate array "+name)
    if(arr.size == 1):
        return float(arr[0]),1.
    diff = np.diff(arr)
    if(not np.allclose(diff,diff[0],rtol=tol_uniform,atol=1.e-12)):
        raise ValueError("non-uniform coordinate array "+name+" not supported by ImageData")
    return float(arr[0]),float(diff[0])
def write_vtk(bpdir,fieldnames,x,y,z):
    #
    # ParaView's VTX reader reverses the ADIOS2/Fortran extent into VTK dimensions
    #
    extent_coordinates = [z,y,x]
    origin_x,spacing_x = infer_origin_spacing(x,"x")
    origin_y,spacing_y = infer_origin_spacing(y,"y")
    origin_z,spacing_z = infer_origin_spacing(z,"z")
    extent = "{} {} {} {} {} {}".format(*sum(([0,max(len(coordinate)-1,0)] for coordinate in extent_coordinates),[]))
    VTKFile = Element("VTKFile",attrib={"type":"ImageData","version":"0.1","byte_order":"LittleEndian"})
    image = SubElement(VTKFile,"ImageData",attrib={"WholeExtent":extent, \
                                                   "Origin":"{:0.16E} {:0.16E} {:0.16E}".format(origin_x,origin_y,origin_z), \
                                                   "Spacing":"{:0.16E} {:0.16E} {:0.16E}".format(spacing_x,spacing_y,spacing_z)})
    piece = SubElement(image,"Piece",attrib={"Extent":extent})
    pointdata = SubElement(piece,"PointData")
    if(len(fieldnames) > 0):
        pointdata.set("Scalars",fieldnames[0])
    for fieldname in fieldnames:
        dataarray = SubElement(pointdata,"DataArray",attrib={"Name":fieldname})
        dataarray.text = ""
    timedata = SubElement(pointdata,"DataArray",attrib={"Name":"TIME"})
    timedata.text = "\n        time\n      "
    output = ElementTree.tostring(VTKFile,'utf-8')
    output = minidom.parseString(output)
    output = output.toprettyxml(indent="    ",newl='\n')
    vtkfile = open(bpdir/"vtk.xml",'w')
    vtkfile.write(output)
    vtkfile.close()
def main():
    logfile = input("Name of the log file written by CaNS [log_visu_3d.out]: ") or "log_visu_3d.out"
    files = get_unique_files(logfile)
    for filename in files:
        bpdir = Path(filename)
        fieldnames,x,y,z = get_bp_data(bpdir)
        write_vtk(bpdir,fieldnames,x,y,z)
if __name__ == "__main__":
    main()
