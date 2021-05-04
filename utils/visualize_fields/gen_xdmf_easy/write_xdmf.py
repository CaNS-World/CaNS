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
dtype_saves = np.dtype([                                                   \
                        ('file' , 'U100'), ('variable', 'U100'),           \
                        ('imin' , int)   , ('jmin' , int), ('kmin' , int), \
                        ('imax' , int)   , ('jmax' , int), ('kmax' , int), \
                        ('istep', int)   , ('jstep', int), ('kstep', int), \
                        ('time', float)  , ('isave', int)                  \
                       ])
geofile  = "geometry.out"
logfile  = input("Name of the log file written by CaNS [log_visu_3d.out]: ") or "log_visu_3d.out"
gridname = input("Name to be appended to the grid files to prevent overwriting []: ") or ""
xgridfile = "x"+gridname+'.bin'
ygridfile = "y"+gridname+'.bin'
zgridfile = "z"+gridname+'.bin'
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
nmin  = np.array([saves['imin' ][0], saves['jmin' ][0], saves['kmin' ][0]])
nmax  = np.array([saves['imax' ][0], saves['jmax' ][0], saves['kmax' ][0]])
nstep = np.array([saves['istep'][0], saves['jstep'][0], saves['kstep'][0]])
n = ((nmax-nmin+1)/nstep).astype(int)
#
# retrieve some computational parameters
#
data = np.loadtxt(geofile, comments = "!", max_rows = 2)
ng = data[0,:].astype('int')
l  = data[1,:]
dl = l/(1.*ng)
#
# generate grid files
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
x[nmin[0]-1:nmax[0]:nstep[0]].astype('float64').tofile(xgridfile)
y[nmin[1]-1:nmax[1]:nstep[1]].astype('float64').tofile(ygridfile)
z[nmin[2]-1:nmax[2]:nstep[2]].astype('float64').tofile(zgridfile)
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
    dataitem.text += "{:15.6E}".format(saves["time"][ii*nflds]) + " "
for ii in range(nsaves):
    grid_fld = SubElement(grid,"Grid", attrib = {"Name": "T{:7}".format(str(saves['isave'][ii*nflds]).zfill(7)), "GridType": "Uniform"})
    topology = SubElement(grid_fld, "Topology", attrib = {"Reference": "/Xdmf/Domain/Topology[1]"})
    geometry = SubElement(grid_fld, "Geometry", attrib = {"Reference": "/Xdmf/Domain/Geometry[1]"})
    for jj in range(nflds):
        index = ii*nflds+jj
        #
        # if vector, skip second and third components from the loop, and write the three files at once
        #
        if( (saves['variable'][index-1+0].endswith('_X') and saves['variable'][index-1+1].endswith('_Y') and saves['variable'][index-1+2].endswith('_Z') and \
             saves['variable'][index-1+0][0:-2] ==           saves['variable'][index-1+1][0:-2] ==           saves['variable'][index-1+2][0:-2]) or \
            (saves['variable'][index-2+0].endswith('_X') and saves['variable'][index-2+1].endswith('_Y') and saves['variable'][index-2+2].endswith('_Z') and \
             saves['variable'][index-2+0][0:-2] ==           saves['variable'][index-2+1][0:-2] ==           saves['variable'][index-2+2][0:-2]) \
          ): continue
        #
        # vector
        #
        if(saves['variable'][index+0].endswith('_X') and saves['variable'][index+1].endswith('_Y') and saves['variable'][index+2].endswith('_Z') and \
           saves['variable'][index+0][0:-2] ==           saves['variable'][index+1][0:-2] ==           saves['variable'][index+2][0:-2]):
           attribute = SubElement(grid_fld, "Attribute", attrib = {"AttributeType": "Vector", "Name": "{}".format(saves['variable'][index][:-2]), "Center": "Node"})
           attribute = SubElement(attribute, "DataItem", attrib = {"Function": "JOIN($0, $1, $2)", "ItemType": "Function", "Dimensions": "{} {} {} {}".format(n[2], n[1], n[0], 3)})
           for q in range(3):
              dataitem = SubElement(attribute, "DataItem", attrib = {"Format": "Binary", "DataType": "Float", "Precision": "{}".format(iprecision), "Endian": "Native", "Seek": "{}".format(iseek), "Dimensions": "{} {} {}".format(n[2], n[1], n[0])})
              dataitem.text = saves['file'][index+q]
        #
        # scalar
        #
        else:
           attribute = SubElement(grid_fld, "Attribute", attrib = {"AttributeType": "Scalar", "Name": "{}".format(saves['variable'][index]), "Center": "Node"})
           dataitem = SubElement(attribute, "DataItem", attrib = {"Format": "Binary", "DataType": "Float", "Precision": "{}".format(iprecision), "Endian": "Native", "Seek": "{}".format(iseek), "Dimensions": "{} {} {}".format(n[2], n[1], n[0])})
           dataitem.text = saves['file'][index]
output = ElementTree.tostring(Xdmf, 'utf-8')
output = minidom.parseString(output)
output = output.toprettyxml(indent="    ",newl='\n')
#
# write visualization file
#
outfile = input("Name of the output file [viewfld_DNS.xmf]: ") or "viewfld_DNS.xmf"
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
