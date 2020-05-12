import numpy as np
import os
#
# define some custom parameters, not defined in the DNS code
#
iseek      = 0            # number of bytes to skip relative to the origin of the binary file (0 for CaNS)
iprecision = 8            # precision of real-valued data
r0 = np.array([0.,0.,0.]) # domain origin
non_uniform_grid = True
#
# define data type and
# read saved data log
#
dtype_saves = np.dtype([                                                         \
                        ('file' , 'U100'), ('variable', 'U100'),                 \
                        ('imin' , np.int), ('jmin' , np.int), ('kmin' , np.int), \
                        ('imax' , np.int), ('jmax' , np.int), ('kmax' , np.int), \
                        ('istep', np.int), ('jstep', np.int), ('kstep', np.int), \
                        ('time', np.float), ('isave', np.int)                    \
                       ])
geofile  = "geometry.out"
logfile  = input("Name of the log file written by CaNS [log_visu_3d.out]: ") or "log_visu_3d.out"
saves = np.loadtxt(logfile, dtype=dtype_saves)
#
# remove duplicates
#
saves = np.unique(saves)
#
# sort elements by increasing isave
#
saves = np.sort(saves, order='isave')
#
# harvest some information from the log file
#
nelements = saves.size
nflds     = 0
isave = saves['isave'][0]
while(isave == saves['isave'][nflds] and nflds < nelements-1): nflds += 1
if(nflds == nelements-1): nflds += 1
nsaves = int(nelements/nflds)
#
# retrieve some computational parameters
#
data = np.loadtxt(geofile, comments = "!", max_rows = 2)
ng = data[0,:].astype('int')
l  = data[1,:]
dl = l/(1.*ng)
#
# write xml file
#
from xml.etree import ElementTree
from xml.dom import minidom
from xml.etree.ElementTree import Element, SubElement, Comment
Xdmf = Element("Xdmf", attrib = {"xmlns:xi": "http://www.w3.org/2001/XInclude", "Version": "2.0"})
domain = SubElement(Xdmf, "Domain")
topology = SubElement(domain,"Topology", attrib = {"name": "TOPO", "TopologyType": "3DRectMesh", "Dimensions" : "{} {} {}".format(2, 2, 2)})
geometry = SubElement(domain,"Geometry", attrib = {"name": "GEO", "GeometryType": "ORIGIN_DXDYDZ"})
dataitem = SubElement(geometry, "DataItem", attrib = {"Format": "XML", "NumberType": "Float", "Dimensions": "{}".format(3)})
dataitem.text = "{:15.6E} {:15.6E} {:15.6E}".format(r0[0],r0[1],r0[2])
dataitem = SubElement(geometry, "DataItem", attrib = {"Format": "XML", "NumberType": "Float", "Dimensions": "{}".format(3)})
dataitem.text = "{:15.6E} {:15.6E} {:15.6E}".format(l[0],l[1],l[2])
grid = SubElement(domain, "Grid", attrib = {"Name": "TimeSeries", "GridType": "Collection",  "CollectionType": "Temporal"})
time = SubElement(grid, "Time", attrib = {"TimeType":"List"})
dataitem = SubElement(time, "DataItem", attrib = {"Format": "XML", "NumberType": "Float", "Dimensions": "{}".format(nsaves)})
dataitem.text = ""
for ii in range(nsaves):
    dataitem.text += "{:15.6E}".format(saves["time"][ii*nflds]) + " "
for ii in range(nsaves):
    grid_fld = SubElement(grid,"Grid", attrib = {"Name": "T{:7}".format(str(saves['isave'][ii*nflds]).zfill(7)), "GridType": "Uniform"})
    topology = SubElement(grid_fld, "Topology", attrib = {"Reference": "/Xdmf/Domain/Topology[1]"})
    geometry = SubElement(grid_fld, "Geometry", attrib = {"Reference": "/Xdmf/Domain/Geometry[1]"})
output = ElementTree.tostring(Xdmf, 'utf-8')
output = minidom.parseString(output)
output = output.toprettyxml(indent="    ",newl='\n')
#
# write visualization file
#
outfile = input("Name of the output file [viewbox_DNS.xmf]: ") or "viewbox_DNS.xmf"
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
