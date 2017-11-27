#!/bin/python
 
#Copyleft Arvind Ravichandran
#Sun Mar 12 00:30:25 CET 2017
#cells.py
#Description:
 
import numpy as np
from scipy.spatial.distance import pdist

class Cells:
    """ data structure for storing cell information"""

    def __init__(self, xi, ori, sim):

       ### assign centre of mass cell positions

        self.xi = xi
        self.ori = ori

        self.vel = {}

        return

    def calculate_ORI(self, sim):
        """ calculate and average the orientational auto correlation"""

        ndelay = int(sim.nsteps/2.0)
        delay = np.zeros((ndelay), dtype=np.float64)
        msd = np.zeros((ndelay), dtype=np.float64)

        for d in range(1, ndelay):
            delay[d] = d
            msd[d]  = np.mean((self.ori[d:].T*self.ori[:-d].T).sum(axis=0).T.mean(axis=0))

        return delay, msd


    def pdistwithpbc(self, sim, t):
        # Split xi into x and y components
        xx = self.xi[t,:,0];    yy = self.xi[t,:,1]

        # find the correct distance matrix for pbc 
        for (xd,ll) in [(xx,sim.lx),(yy,sim.ly)]:
            # find the periodic 1d distance
            pd = pdist(xd.reshape(xd.shape[0],1))

            # apply boundary condition
            pd = pd-ll*np.round(pd/ll)

            try:
                # sum                                           (if xx)
                total+=pd**2
            except(NameError):
                # or define the sum if not previously defined   (if y)
                total=pd**2

        pd = np.sqrt(total)
        return pd
    
    def pori(self, sim, t):
        # prepare indices of upper triangular matrix without diagonal
        iu = np.triu_indices(sim.nfils,1)    
        orit = self.ori[t,:]
        return dotproduct(orit[iu[0]], orit[iu[1]])

    def dispcorr(self, sim, t, lag):
        # prepare indices of upper triangular matrix without diagonal
        iu = np.triu_indices(sim.nfils,1)    

        # Displacements of all centers of mass
        dx = self.xi[t+lag,:,:] - self.xi[t,:,:]

        # returns dot(di,dj), normalisation (see wysocki1402)
        return dotproduct(dx[iu[0]], dx[iu[1]]), sum(dotproduct(dx,dx))

    def displacements(self, sim, lag):
        dx = self.xi[:-lag, :, :] - self.xi[lag:, :, :]
        self.vel[lag] = dx

    def velocity_projection_not_normalised(self, lag):
        v   = np.array([xx.T for xx in self.vel[lag]])
        ori = self.ori[:-lag] 
        return np.sum(v*ori, axis=2) 

    def velocity_projection_normalised(self, lag):
        v   = np.array([xx.T for xx in self.vel[lag]])
        ori = self.ori[:-lag] 

        normed = np.linalg.norm(v, axis=2)[:,:,np.newaxis]
        vnorm = v/normed

        return np.sum(vnorm*ori, axis=2) 
