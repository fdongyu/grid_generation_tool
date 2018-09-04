# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 14:23:09 2018

@author: dongyu
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 13:27:31 2018

@author: dongyu
"""

import numpy as np
from matplotlib.collections import PolyCollection, LineCollection
import utm
import matplotlib.pyplot as plt


from sunpy import edgeplot


import pdb


class SUNGrid_tools(object):
    """
    general class for ploting a grid given the grid directory
    """    
    _FillValue=999999  

    #originx = 329500
    originx = 616832
    originy = 2854748      
    
    def __init__(self,dire, gridname, **kwargs):
        self.__dict__.update(kwargs)
        
        self.dire = dire
        self.gname = gridname
        
        self.read_grid()
        self.calc_dg()  
        self.dg_plot()
        
        self.grid_plot()
    
    def read_grid(self):
        
        edgedata = self.readTXT('%s/edges.dat'%self.dire)
        self.grad = np.asarray(edgedata[:,3:5],int)
        self.edges = np.asarray(edgedata[:,0:2],int)
        
        celldata = self.readTXT('%s/cells.dat'%self.dire)
        self.cells = np.asarray(celldata[:,2:5],int)
        self.Nc = celldata.shape[0]
        self.nfaces = 3*np.ones((self.Nc,),np.int)
        self.xv = (np.asarray(celldata[:,0], float) - self.originx)/1.
        self.yv = (np.asarray(celldata[:,1], float) - self.originy)/1.
        
        pointdata = self.readTXT('%s/points.dat'%self.dire)
        self.xp = (pointdata[:,0] - self.originx)/1.
        self.yp = (pointdata[:,1] - self.originy)/1.
        
        self.cellmask = self.cells==int(self._FillValue)
        
        self.maxfaces = 3 
        self.Np = len(self.xp)
        

    def calc_dg(self):
        """
        Manually calculate the distance between voronoi points, 'dg'
        """
        print 'Calculating dg...'
        print np.shape(self.grad)
        
        grad = self.grad
        Ne = len(grad)
        for ii in range(Ne):
            if grad[ii,0]==-1:
                grad[ii,0]=grad[ii,1]
            elif grad[ii,1]==-1:
                grad[ii,1]=grad[ii,0]
                
                
        x1 = self.xv[grad[:,0]]
        x2 = self.xv[grad[:,1]]
        y1 = self.yv[grad[:,0]]
        y2 = self.yv[grad[:,1]]
        
        dx=x1-x2
        dy=y1-y2
        
        self.dg = np.sqrt( dx*dx + dy*dy )
        
        #pdb.set_trace()
        

    def badedges(self, threshold):
        """
        find the index of bad edges
        """
        edgedata = self.readTXT('%s/edges.dat'%self.dire)
        
        edge = np.asarray(edgedata[:,0:2], np.int)
        marker = np.asarray(edgedata[:,2], np.int)
        
        ind_boundary = []
        for i in range(len(marker)):
            if marker[i] != 0:
                ind_boundary.append(i)
        
        ind_tem = np.squeeze(np.argwhere(self.dg<threshold)).tolist()
        ind = list(set(ind_tem)-set(ind_boundary))
        
        return np.asarray(ind)
        
        
    
    def dg_plot(self):
        
        xylines = [self.xp[self.edges],self.yp[self.edges]]        
        
        
        fig = plt.figure(figsize=(10,15))
        ax = fig.add_subplot(111)
        
        clim = [0,500]
        xlim = np.asarray([616832, 709049]) - self.originx
        ylim = np.asarray([2854748, 3050710]) - self.originy    
    
        # Find the colorbar limits if unspecified
        if clim==None:
            clim=[]
            clim.append(np.min(self.dg))
            clim.append(np.max(self.dg))
     
        # Create the inputs needed by line collection
        Ne = xylines[0].shape[0]
    
        # Put this into the format needed by LineCollection
        linesc = [zip(xylines[0][ii,:],xylines[1][ii,:]) for ii in range(Ne)]    
        collection = LineCollection(linesc,array=self.dg, cmap='bone', linewidth=0.2)
        collection.set_clim(vmin=clim[0],vmax=clim[1])        
        ax.add_collection(collection)
        
        ## extremely small dg
        ind = self.badedges(threshold=50)
        print "Bad edge indexes are ...\n"
        print ind
        #ind = np.squeeze(np.argwhere(self.dg<50))
        #pdb.set_trace()
        linesc2 = [zip(xylines[0][ii,:],xylines[1][ii,:]) for ii in ind]  
        collection2 = LineCollection(linesc2,colors='r', linewidths=2.5)
        ax.add_collection(collection2)
        
        ## plot outline
        #lc = self.outline() #line collection of outline
        #ax.add_collection(lc)
        
        ax.set_aspect('equal')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        cb = fig.colorbar(collection, orientation='vertical')        
        
        plt.show()
        
        
        
        
    def grid_plot(self):
        """
        plotting the grid
        """            
        maxfaces = self.cells.shape[1]
        
        pointdata = self.readTXT('%s/points.dat'%self.dire)
        xp1 = (pointdata[:,0] - self.originx)/1.
        yp1 = (pointdata[:,1] - self.originy)/1.
        #########################################################
        xp = np.zeros((self.Nc,maxfaces+1))
        yp = np.zeros((self.Nc,maxfaces+1))
            
        cells=self.cells.copy()
        #cells[cells.mask]=0
        xp[:,:maxfaces]=xp1[cells]
        xp[range(self.Nc),self.nfaces]=xp1[cells[:,0]]
        yp[:,:maxfaces]=yp1[cells]
        yp[range(self.Nc),self.nfaces]=yp1[cells[:,0]]
        ##########################################################
        xy = np.zeros((maxfaces+1,2))
        def _closepoly(ii):
            nf=self.nfaces[ii]+1
            xy[:nf,0]=xp[ii,:nf]
            xy[:nf,1]=yp[ii,:nf]
            return xy[:nf,:].copy()

        cellxy= [_closepoly(ii) for ii in range(self.Nc)]
        
        
        #xlims = np.asarray([self.xv.min(), self.xv.max()]) 
        #ylims = np.asarray([self.yv.min(), self.yv.max()])   
        
        xlims = np.asarray([616832, 709049]) - self.originx
        ylims = np.asarray([2854748, 3050710]) - self.originy  
        
        
        fig = plt.figure(figsize=(12,10))
        ax = fig.add_subplot(111)
        #cmap = plt.set_cmap('bwr')
        collection = PolyCollection(cellxy,facecolors='none')
        #collection.set_array(np.array(self.L[:]))
        collection.set_edgecolors('k')
        collection.set_linewidths(0.3)
        ax.add_collection(collection)
        lc = self.outline() #line collection of outline
        ax.add_collection(lc)
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        ax.set_aspect('equal')
        ax.tick_params(labelsize=22)
        ax.set_xlabel('Easting (m)', fontsize=22)
        ax.set_ylabel('Northing (m)', fontsize=22)
        fig.tight_layout()              
    
        plt.show()
        #plt.savefig('figures/%s.png'%self.gname, format='png')
        #plt.show()
        #plt.close()   
        
        
        
    def outline(self):
        """
        function that used to plot the outline of the grid
        """
        edgedata = self.readTXT('%s/edges.dat'%self.dire)
        
        edge = np.asarray(edgedata[:,0:2], np.int)
        marker = np.asarray(edgedata[:,2], np.int)
        
        edge_new = []
        for i in range(len(marker)):
            if marker[i] != 0:
                edge_new.append(edge[i])
        
        edge_new = np.asarray(edge_new)

        pointdata = self.readTXT('%s/points.dat'%self.dire)
        xp = (pointdata[:,0] - self.originx)/1.
        yp = (pointdata[:,1] - self.originy)/1.
        
        lines = []
        for j in range(edge_new.shape[0]):
            ind1 = edge_new[j,0]
            ind2 = edge_new[j,1]
            lines.append([(xp[ind1],yp[ind1]),(xp[ind2],yp[ind2])])

        return LineCollection(lines,colors='k', linewidths=1.8) 
        
    def readTXT(self,fname,sep=None):
        """
        Reads a txt file into an array of floats
        """
        
        fp = file(fname,'rt')
        data = np.array([map(float,line.split(sep)) for line in fp])
        fp.close()
        
        return data
 


       
#### For testing ####
if __name__ == "__main__":

    dire = '../Laguna_Madre/grid_generation/coarse_new2/grids'
    gname = 'LM_coarse'
    SUNGrid_tools(dire, gname)



