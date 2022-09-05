#!/usr/bin/env python
#
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#
def initfig(width = 600.0, ar = (np.sqrt(5)-1.0)/2.0, scl = 0.95):
    fig_width_pt  = width                                # \showthe\columnwidth
    inches_per_pt = 1.0/72.27                            # pt to in
    aspect_ratio  = ar                                   # aspect ratio
    fig_scale     = scl                                  # scale factor
    fig_width     = fig_width_pt*inches_per_pt*fig_scale # width in in
    fig_height    = fig_width*aspect_ratio               # height in in
    fig_size      = [fig_width,fig_height]               # final dimensions
    params = {'backend'        : 'ps',
              'font.family'    : 'serif',
              'font.size'      : 10,
              'axes.labelsize' : 10,
              'legend.fontsize': 8,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'text.usetex'    : True,
              'figure.figsize' : fig_size}
    plt.rcParams.update(params)
#
def format(x, pos):
    'The two args are the value and tick position'
    return '%1.2g' % (x)
#
def newfig():
    plt.cla()
    plt.clf()
    plt.close()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    formatter = matplotlib.ticker.FuncFormatter(format)
    return fig, ax, formatter
