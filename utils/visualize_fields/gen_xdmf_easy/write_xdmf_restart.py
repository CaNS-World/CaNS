#
# this python script can be used to visualize directly field from the checkpoint files of CaNS
#
import numpy as np
import os
import glob
import struct
#
# define some custom parameters, not defined in the DNS code
#
iprecision = 8            # precision of the real-valued data
r0 = np.array([0.,0.,0.]) # domain origin
non_uniform_grid = True
#
# retrieve input information
#
filenames = input("Name of the pattern of the restart files to be visualzied [fld?*.bin]: ") or "fld?*.bin"
files     = glob.glob(filenames)
nsaves    = np.size(files)
variables = input("Names of stored variables [VEX VEY VEZ PRE]: ") or "VEX VEY VEZ PRE"
variables = variables.split(" ")
nflds     = np.size(variables)
gridname  = input("Name to be appended to the grid files to prevent overwriting [_fld]: ") or "_fld"
xgridfile = "x"+gridname+'.bin'
ygridfile = "y"+gridname+'.bin'
zgridfile = "z"+gridname+'.bin'
#
# retrieve other computational parameters
#
geofile   = "geometry.out"
data = np.loadtxt(geofile, comments = "!", max_rows = 2)
ng = data[0,:].astype('int')
l  = data[1,:]
dl = l/(1.*ng)
n  = ng
rtimes = np.zeros(nsaves)
isteps = np.zeros(nsaves,dtype=int)
iseek   = n[0]*n[1]*n[2]*iprecision*nflds # file offset in bytes with respect to the origin
                                          # (to retrieve the simulation time and time step number)
for i in range(nsaves):
    with open(files[i], 'rb') as f:
        raw   = f.read()[iseek:iseek+iprecision*2]
    rtimes[i] =     struct.unpack('2d',raw)[0]
    isteps[i] = int(struct.unpack('2d',raw)[1])
    f.close()
#
# remove duplicates
#
isteps, indeces = np.unique(isteps,return_index=True)
rtimes  = np.take(rtimes, indeces)
files   = np.take(files,  indeces)
nsaves  = np.size(files)
#
# sort by increasing istep
#
indeces = np.argsort(isteps)
isteps  = np.take(isteps, indeces)
rtimes  = np.take(rtimes, indeces)
files   = np.take(files,  indeces)
#
# create grid files
#
x = np.arange(r0[0]+dl[0]/2.,r0[0]+l[0],dl[0])
y = np.arange(r0[1]+dl[1]/2.,r0[1]+l[1],dl[1])
z = np.arange(r0[2]+dl[2]/2.,r0[2]+l[2],dl[2])
if os.path.exists(xgridfile): os.remove(xgridfile)
if os.path.exists(ygridfile): os.remove(ygridfile)
if os.path.exists(zgridfile): os.remove(zgridfile)
if(non_uniform_grid):
    f   = open('grid.bin','rb')
    if(    iprecision == 4):
        grid_z = np.fromfile(f,dtype='float32')
    else:
        grid_z = np.fromfile(f,dtype='float64')
    f.close()
    grid_z = np.reshape(grid_z,(ng[2],4),order='F')
    z = r0[2] + grid_z[:,2]
x[0:n[0]].astype('float64').tofile(xgridfile)
y[0:n[1]].astype('float64').tofile(ygridfile)
z[0:n[2]].astype('float64').tofile(zgridfile)
#
# write xml file
#
from xml.etree import ElementTree
from xml.dom import minidom
from xml.etree.ElementTree import Element, SubElement, Comment
Xdmf = Element("Xdmf", attrib = {"xmlns:xi": "http://www.w3.org/2001/XInclude", "Version": "2.0"})
domain = SubElement(Xdmf, "Domain")
topology = SubElement(domain,"Topology", attrib = {"name": "TOPO", "TopologyType": "3DRectMesh", "Dimensions" : "{} {} {}".format(n[2], n[1], n[0])})
geometry = SubElement(domain,"Geometry", attrib = {"name": "GEO", "GeometryType": "VXVYVZ"})
dataitem = SubElement(geometry, "DataItem", attrib = {"Format": "Binary", "DataType": "Float", "Precision": "{}".format(iprecision), "Endian": "Native", "Dimensions": "{}".format(n[0])})
dataitem.text = xgridfile
dataitem = SubElement(geometry, "DataItem", attrib = {"Format": "Binary", "DataType": "Float", "Precision": "{}".format(iprecision), "Endian": "Native", "Dimensions": "{}".format(n[1])})
dataitem.text = ygridfile
dataitem = SubElement(geometry, "DataItem", attrib = {"Format": "Binary", "DataType": "Float", "Precision": "{}".format(iprecision), "Endian": "Native", "Dimensions": "{}".format(n[2])})
dataitem.text = zgridfile
grid = SubElement(domain, "Grid", attrib = {"Name": "TimeSeries", "GridType": "Collection",  "CollectionType": "Temporal"})
time = SubElement(grid, "Time", attrib = {"TimeType":"List"})
dataitem = SubElement(time, "DataItem", attrib = {"Format": "XML", "NumberType": "Float", "Dimensions": "{}".format(nsaves)})
dataitem.text = ""
for ii in range(nsaves):
    dataitem.text += "{:15.6E}".format(rtimes[ii]) + " "
for ii in range(nsaves):
    grid_fld = SubElement(grid,"Grid", attrib = {"Name": "T{:7}".format(str(isteps[ii]).zfill(7)), "GridType": "Uniform"})
    topology = SubElement(grid_fld, "Topology", attrib = {"Reference": "/Xdmf/Domain/Topology[1]"})
    geometry = SubElement(grid_fld, "Geometry", attrib = {"Reference": "/Xdmf/Domain/Geometry[1]"})
    for jj in range(nflds):
        iseek = jj*iprecision*n[2]*n[1]*n[0]
        attribute = SubElement(grid_fld, "Attribute", attrib = {"Name": "{}".format(variables[jj]), "Center": "Node"})
        dataitem = SubElement(attribute, "DataItem", attrib = {"Format": "Binary", "DataType": "Float", "Precision": "{}".format(iprecision), "Endian": "Native", "Seek": "{}".format(iseek), "Dimensions": "{} {} {}".format(n[2], n[1], n[0])})
        dataitem.text = files[ii]
output = ElementTree.tostring(Xdmf, 'utf-8')
output = minidom.parseString(output)
output = output.toprettyxml(indent="    ",newl='\n')
#
# write visualization file
#
outfile = input("Name of the output file [viewfld_DNS_fld.xmf]: ") or "viewfld_DNS_fld.xmf"
xdmf_file = open(outfile, 'w')
xdmf_file.write(output)
xdmf_file.close
#
# workaround to add the DOCTYPE line
#
xdmf_file = open(outfile, "r")
contents = xdmf_file.readlines()
xdmf_file.close()
contents.insert(1, '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
xdmf_file = open(outfile, "w")
contents = "".join(contents)
xdmf_file.write(contents)
xdmf_file.close()
