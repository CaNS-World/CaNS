#!/usr/bin/env python
from pltparams import *
from param import *
import pylab
#
# read data
#
f   = open(filenamei,'rb')
fld = np.fromfile(f,dtype='float64')
fld = np.reshape(fld,(n2,n1),order='C')
f.close()
#
# initialize figure
#
initfig(ar = l2/l1)
fig, ax, formatter = newfig()
#ax.set_axis_off()
#ax.get_xaxis().set_visible(False)
#ax.get_yaxis().set_visible(False)
#
# plot data
#
x1 = np.linspace(0.+dx1/2.,l1-dx1/2.,n1)
x2 = np.linspace(0.+dx2/2.,l2-dx2/2.,n2)
cs1 = ax.contourf(x1, x2, fld,
                  cmap=plt.cm.jet)
cbar = fig.colorbar(cs1, orientation='vertical')
#
# format figrue
#
ax.set_title(fgtitle)
cbar.ax.set_title(cbtitle)
ax.set_xlabel(x1title)
ax.set_ylabel(x2title)
#
# save figure
#
fig.tight_layout(pad=0.15)
plt.show()
fig.savefig(filenameo)
