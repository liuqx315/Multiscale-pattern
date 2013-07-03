#!/usr/bin/env python
# matplotlib-based plotting script for heat1D.c example
# Daniel R. Reynolds, reynolds@smu.edu

# imports
import sys
import pylab as plt
import numpy as np

# load mesh data file
inp = open('heat_mesh.txt').readlines()
mesh = []
for line in inp:
    mesh.append(np.array(str.split(line), dtype=np.double))


# load solution data file
inp = open('heat1D.txt').readlines()
data = []
for line in inp:
    data.append(np.array(str.split(line), dtype=np.double))

# determine number of time steps
nt  = len(mesh)
nt2 = len(data)
if (nt != nt2):
    sys.exit('plot_heat1D_adapt.py error: data and mesh files have different numbers of time steps')

# determine maximum temperature
maxtemp = 0.0
for tstep in range(nt):
    mx = 1.1*data[tstep].max()
    if (mx > maxtemp):
        maxtemp = mx

# generate plots of results
for tstep in range(nt):

    # set string constants for output plots, current time, mesh size
    pname = 'heat1d.' + repr(tstep).zfill(3) + '.png'
    tstr  = repr(tstep)
    nxstr = repr(len(data[tstep]))

    # plot current solution and save to disk
    plt.figure(1)
    plt.plot(mesh[tstep],data[tstep])
    plt.xlabel('x')
    plt.ylabel('solution')
    plt.title('u(x) at output ' + tstr + ', mesh = ' + nxstr)
    plt.axis((0.0, 1.0, 0.0, maxtemp))
    plt.grid()
    plt.savefig(pname)
    plt.close()


##### end of script #####
