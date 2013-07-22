#!/usr/bin/env python
# matplotlib-based plotting script for heat2D.cpp example
# Daniel R. Reynolds, reynolds@smu.edu

# imports
import sys
import numpy as np
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

rou_data = np.loadtxt('2DWENO_rou.txt', dtype=np.double)
qx_data = np.loadtxt('2DWENO_qx.txt', dtype=np.double)
qy_data = np.loadtxt('2DWENO_qy.txt', dtype=np.double)
E_data = np.loadtxt('2DWENO_E.txt', dtype=np.double)
nt = np.shape(rou_data)[0];

xspan = np.loadtxt('WENO2D_meshx.txt', dtype=np.double)
yspan = np.loadtxt('WENO2D_meshy.txt', dtype=np.double)

ny = 400
nx = 400

results_rou = np.zeros((nt,ny,nx));
results_qx = np.zeros((nt,ny,nx));
results_qy = np.zeros((nt,ny,nx));
results_E = np.zeros((nt,ny,nx));

for i in range(nt):
    results_rou[i,:,:] = np.reshape(rou_data[i,:], (ny,nx))
for i in range(nt):
    results_qx[i,:,:] = np.reshape(qx_data[i,:], (ny,nx))
for i in range(nt):
    results_qy[i,:,:] = np.reshape(qy_data[i,:], (ny,nx))
for i in range(nt):
    results_E[i,:,:] = np.reshape(E_data[i,:], (ny,nx))

# determine extents of plots
maxtemp_rou = 1.1*results_rou.max()
mintemp_rou = 0.9*results_rou.min()
maxtemp_qx = 1.1*results_qx.max()
mintemp_qx = 0.9*results_qx.min()
maxtemp_qy = 1.1*results_qy.max()
mintemp_qy = 0.9*results_qy.min()
maxtemp_E = 1.1*results_E.max()
mintemp_E = 0.9*results_E.min()

maxtemp = np.array([maxtemp_rou, maxtemp_qx, maxtemp_qy, maxtemp_E]).max();
mintemp = np.array([mintemp_rou, mintemp_qx, mintemp_qy, mintemp_E]).min();

for tstep in range(nt):

# set string constants for output plots, current time, mesh size
    pname = 'WENO2d_surf.' + repr(tstep).zfill(3) + '.png'
    cname = 'WENO2d_contour.' + repr(tstep).zfill(3) + '.png'
    tstr  = repr(tstep)
    nxstr = repr(nx)
    nystr = repr(ny)

X,Y = np.meshgrid(xspan,yspan)

# plot current solution as a surface, and save to disk
    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, results_rou[tstep,:,:], rstride=1, cstride=1, 
                    cmap=cm.jet, linewidth=0, antialiased=True, shade=True)
    ax.plot_surface(X, Y, results_qx[tstep,:,:], rstride=1, cstride=1, 
                    cmap=cm.jet, linewidth=0, antialiased=True, shade=True)
    ax.plot_surface(X, Y, results_qy[tstep,:,:], rstride=1, cstride=1, 
                    cmap=cm.jet, linewidth=0, antialiased=True, shade=True)
    ax.plot_surface(X, Y, results_E[tstep,:,:], rstride=1, cstride=1, 
                    cmap=cm.jet, linewidth=0, antialiased=True, shade=True)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlim((mintemp, maxtemp))
    ax.view_init(20,45)
    title('u(x,y) at output ' + tstr + ', mesh = ' + nxstr + 'x' + nystr)
    savefig(pname)
    plt.close()


    plt.contourf(xspan,yspan,results_rou[tstep,:,:],15, cmap=plt.cm.jet)
    plt.contourf(xspan,yspan,results_qx[tstep,:,:],15, cmap=plt.cm.jet)
    plt.contourf(xspan,yspan,results_qy[tstep,:,:],15, cmap=plt.cm.jet)
    plt.contourf(xspan,yspan,results_E[tstep,:,:],15, cmap=plt.cm.jet)
    plt.colorbar()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('contour at output ' + tstr + ', mesh = ' + nxstr + 'x' + nystr)
    plt.savefig(cname)
    plt.close()
