# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 11:02:54 2016

@author: dongyu
"""

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
import os
import string

import pdb

fname = os.getcwd() + '/test/Bath-Laguna.2dm'
#fname = os.getcwd() + '/test/Bath-Matag.2dm'
fname = 'ALL_BAYS_latlon_NAVD88_meters.2dm'

## Read txt file
data = []
a = open(fname, 'r').readlines()
for s in a:
    line = s.split()
    data.append(line)

data1 = []
data2 = []
for l in data:
    if len(l) == 6:
        data1.append(l)
    elif len(l) == 5:
        data2.append(l)
    else:
        print "No such option!!!"

## initialize the matrix 
cells = np.zeros([len(data1), 3])
points = np.zeros([len(data2), 2])
depth = np.zeros([len(data2)])

for i in range(len(data1)):
    cells[i,:] = [string.atoi(data1[i][2]), string.atoi(data1[i][3]), string.atoi(data1[i][4])]
    
for j in range(len(data2)):
    points[j,:] = [string.atof(data2[j][2]), string.atof(data2[j][3])]
    depth[j] = string.atof(data2[j][4])
    
## collection plotting
Nc = len(data1)
nfaces = 3*np.ones((Nc,),np.int)
MAXFACES = cells.shape[1]
xp1 = points[:,0]
yp1 = points[:,1]

maxfaces=MAXFACES
xp = np.zeros((Nc,maxfaces+1))
yp = np.zeros((Nc,maxfaces+1))

## important for sms type of grid: start from 1
cells = cells - 1
cells = np.asarray(cells, int)
cells = cells.copy()
#pdb.set_trace()
xp[:,:maxfaces]=xp1[cells]
xp[range(Nc),nfaces]=xp1[cells[:,0]]
yp[:,:maxfaces]=yp1[cells]
yp[range(Nc),nfaces]=yp1[cells[:,0]]

xlims=(xp.min(),xp.max())
ylims=(yp.min(),yp.max())

xy = np.zeros((maxfaces+1,2))
def _closepoly(ii):
    nf=nfaces[ii]+1
    xy[:nf,0]=xp[ii,:nf]
    xy[:nf,1]=yp[ii,:nf]
    return xy[:nf,:].copy()

cellxy= [_closepoly(ii) for ii in range(Nc)]

dep=[]
for i in range(Nc):
    dep.append((depth[cells[i][0]]+depth[cells[i][1]]+depth[cells[i][2]])/3)
    
clim=[min(dep),max(dep)]
fig=plt.figure(figsize=(15,15))
axes = fig.add_subplot(111)
collection = PolyCollection(cellxy,facecolors='none',cmap='jet')
#collection.set_array(np.array(dep[:]))
collection.set_edgecolors('k')
collection.set_linewidths(0.2)
#collection.set_clim(vmin=clim[0],vmax=clim[1])
#collection.set_edgecolors(collection.to_rgba(np.array((dep[:])))) 
axes.add_collection(collection)
axes.set_xlim(xlims)
axes.set_ylim(ylims)
axes.set_aspect('equal')
axes.set_xlabel('Easting [m]')
axes.set_ylabel('Northing [m]')
plt.title('Laguna Madre Bay')
#plt.show() 
plt.savefig('All_Texas_Bays.png')
plt.close()


pdb.set_trace()
