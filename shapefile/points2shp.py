# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 10:55:47 2018

@author: dongyu
"""

import numpy as np
import osgeo.ogr, osgeo.osr
import os
import string
import utm

import pdb



#### read data
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
    points[j,:] = [string.atof(data2[j][3]), string.atof(data2[j][2])]
    depth[j] = string.atof(data2[j][4])
      
xp = points[:,0]
yp = points[:,1]

sw = [25.648536, -98.088590]
ne = [27.933344, -96.638395]

ll = np.array([sw[0],sw[1]])
ur = np.array([ne[0],ne[1]])

inidx = np.all(np.logical_and(ll <= points, points<=ur), axis=1)

points_new = points[inidx]

#### create the shape file
shpfile = 'Laguna_Madre.shp'
srs = osgeo.osr.SpatialReference()
srs.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
#proj = "UTM %d (%s) in northern hemisphere."%(zone,CS)

driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
if os.path.exists(shpfile):
    os.unlink(shpfile)
shapeData = driver.CreateDataSource(shpfile)
layer = shapeData.CreateLayer('Grid', srs, osgeo.ogr.wkbPoint)
layerDefinition = layer.GetLayerDefn()
# Loop through the list of nodes to get the coordinates of each polygon
#multipoint = osgeo.ogr.Geometry(osgeo.ogr.wkbMultiPoint)
ctr=0
for xy in points_new:
    ctr+=1
    point = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
    point.SetPoint(0, xy[1], xy[0])
    
    featureIndex = ctr
    feature = osgeo.ogr.Feature(layerDefinition)
    feature.SetGeometry(point)
    feature.SetFID(featureIndex)
    layer.CreateFeature(feature)
    feature.Destroy()    

#ctr=0
#for xy in points_new:
#    ctr+=1
#    
#    # Add points individually to the polygon
#    #for nodes in xy:
#    
#    ring.AddPoint(xy[0],xy[1])
#            
#    poly = osgeo.ogr.Geometry(osgeo.ogr.wkbPolygon)
#    poly.AddGeometry(ring)
#    
#    # Update the feature with the polygon data
#    featureIndex = ctr
#    feature = osgeo.ogr.Feature(layerDefinition)
#    feature.SetGeometry(poly)
#    feature.SetFID(featureIndex)
#    layer.CreateFeature(feature)
#    feature.Destroy()        

shapeData.Destroy()
print 'Complete - file written to:\n      %s'%shpfile

pdb.set_trace()











