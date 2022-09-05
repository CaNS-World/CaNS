#!/usr/bin/env python
#
datadir = 'data/'
figsdir = 'figs/'
exto = '.pdf'
filenamei = datadir + 'fld_u_slice_fld_0629500.bin'
filenameo = figsdir + 'visu' + exto
#
n2 = 128
n1 = 128
l2 = 1.0
l1 = 1.0
dx1 = l1/float(n1)
dx2 = l2/float(n2)
#
fgtitle = r'Slice at $y/h = 0.5$'
cbtitle = r'$u/U_b$'
x1title = r'$x/h$'
x2title = r'$z/h$'
