# -
#
# SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
# SPDX-License-Identifier: MIT
#
# -
#
# this python script can be used to visualize directly fields from HDF5 checkpoint files of CaNS
#
import glob
import os
import h5py
import numpy as np
#
# define some custom parameters, not defined in the DNS code
#
iprecision = 8
if(iprecision == 4):
    my_dtype = 'float32'
else:
    my_dtype = 'float64'
r0 = np.array([0.,0.,0.]) # domain origin
#
# retrieve input information
#
filenames = input("Name of the pattern of the HDF5 restart files to be visualized [fld?*.h5]: ") or "fld?*.h5"
files = glob.glob(filenames)
nsaves = np.size(files)
variables = input("Names of displayed variables [VEX VEY VEZ PRE]: ") or "VEX VEY VEZ PRE"
variables = variables.split(" ")
fieldnames = ['u','v','w','p']
gridname  = input("Name to be appended to the grid files to prevent overwriting [_fld]: ") or "_fld"
gridfile = "grid"+gridname+'.h5'
#
# retrieve other computational parameters
#
geofile   = "geometry.out"
data = np.loadtxt(geofile, comments = "!", max_rows = 2)
ng = data[0,:].astype('int')
l  = data[1,:]
dl = l/(1.*ng)
n  = ng
def get_split_groups(files):
    groups = {}
    for filename in files:
        root, ext = os.path.splitext(filename)
        if(ext != ".h5"): continue
        for fieldname in fieldnames:
            suffix = "_" + fieldname
            if(root.endswith(suffix)):
                groups.setdefault(root[:-len(suffix)], {})[fieldname] = filename
    return [(prefix, groups[prefix]) for prefix in sorted(groups) \
            if all(fieldname in groups[prefix] for fieldname in fieldnames)]
split_groups = get_split_groups(files)
is_split = len(split_groups) > 0
if(is_split):
    nsaves = len(split_groups)
rtimes = np.zeros(nsaves)
isteps = np.zeros(nsaves,dtype=int)
for i in range(nsaves):
    filename = split_groups[i][1]['u'] if is_split else files[i]
    hf = h5py.File(filename,'r')
    rtimes[i] = float(np.asarray(hf["meta/time"])[0])
    isteps[i] = int(np.asarray(hf["meta/istep"])[0])
    hf.close()
#
# remove duplicates
#
isteps, indeces = np.unique(isteps,return_index=True)
rtimes  = np.take(rtimes, indeces)
if(is_split):
    split_groups = [split_groups[i] for i in indeces]
else:
    files = np.take(files, indeces)
nsaves  = np.size(files)
if(is_split): nsaves = len(split_groups)
#
# sort by increasing istep
#
indeces = np.argsort(isteps)
isteps  = np.take(isteps, indeces)
rtimes  = np.take(rtimes, indeces)
if(is_split):
    split_groups = [split_groups[i] for i in indeces]
else:
    files = np.take(files, indeces)
#
# create grid files
#
x = np.linspace(r0[0]+dl[0]/2.,r0[0]+l[0]-dl[0]/2.,ng[0])
y = np.linspace(r0[1]+dl[1]/2.,r0[1]+l[1]-dl[1]/2.,ng[1])
z = np.linspace(r0[2]+dl[2]/2.,r0[2]+l[2]-dl[2]/2.,ng[2])
if(os.path.exists('grid.h5')):
    hf = h5py.File('grid.h5','r')
    z = np.asarray(hf['rc'])
    hf.close()
elif(os.path.exists('grid.bin')):
    f   = open('grid.bin','rb')
    grid_z = np.fromfile(f,dtype=my_dtype)
    f.close()
    grid_z = np.reshape(grid_z,(ng[2],4),order='F')
    z = r0[2] + grid_z[:,2]
if os.path.exists(gridfile): os.remove(gridfile)
hf = h5py.File(gridfile,'w')
hf.create_dataset('x',data=x.astype(my_dtype))
hf.create_dataset('y',data=y.astype(my_dtype))
hf.create_dataset('z',data=z.astype(my_dtype))
hf.close()
#
# write xml file
#
from xml.etree import ElementTree
from xml.dom import minidom
from xml.etree.ElementTree import Element, SubElement
Xdmf = Element("Xdmf", attrib = {"xmlns:xi": "http://www.w3.org/2001/XInclude", "Version": "2.0"})
domain = SubElement(Xdmf, "Domain")
grid = SubElement(domain, "Grid", attrib = {"Name": "TimeSeries", "GridType": "Collection",  "CollectionType": "Temporal"})
for ii in range(nsaves):
    grid_fld = SubElement(grid,"Grid", attrib = {"Name": "T{:7}".format(str(isteps[ii]).zfill(7)), "GridType": "Uniform"})
    time = SubElement(grid_fld, "Time", attrib = {"Value":"{:15.6E}".format(rtimes[ii])})
    topology = SubElement(grid_fld,"Topology", attrib = {"TopologyType": "3DRectMesh", "Dimensions" : "{} {} {}".format(n[2], n[1], n[0])})
    geometry = SubElement(grid_fld,"Geometry", attrib = {"GeometryType": "VXVYVZ"})
    dataitem = SubElement(geometry, "DataItem", attrib = {"Format": "HDF", "NumberType": "Float", "Precision": "{}".format(iprecision), "Dimensions": "{}".format(n[0])})
    dataitem.text = gridfile + ':/x'
    dataitem = SubElement(geometry, "DataItem", attrib = {"Format": "HDF", "NumberType": "Float", "Precision": "{}".format(iprecision), "Dimensions": "{}".format(n[1])})
    dataitem.text = gridfile + ':/y'
    dataitem = SubElement(geometry, "DataItem", attrib = {"Format": "HDF", "NumberType": "Float", "Precision": "{}".format(iprecision), "Dimensions": "{}".format(n[2])})
    dataitem.text = gridfile + ':/z'
    for jj in range(min(len(variables),len(fieldnames))):
        filename = split_groups[ii][1][fieldnames[jj]] if is_split else files[ii]
        attribute = SubElement(grid_fld, "Attribute", attrib = {"Name": "{}".format(variables[jj]), "Center": "Node"})
        dataitem = SubElement(attribute, "DataItem", attrib = {"Format": "HDF", "NumberType": "Float", "Precision": "{}".format(iprecision), "Dimensions": "{} {} {}".format(n[2], n[1], n[0])})
        dataitem.text = filename + ':/fields/' + fieldnames[jj]
output = ElementTree.tostring(Xdmf, 'utf-8')
output = minidom.parseString(output)
output = output.toprettyxml(indent="    ",newl='\n')
#
# write visualization file
#
outfile = input("Name of the output file [viewfld_DNS.xmf]: ") or "viewfld_DNS.xmf"
xdmf_file = open(outfile, 'w')
xdmf_file.write(output)
xdmf_file.close()
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
