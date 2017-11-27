#!/bin/python
 
#Copyleft Arvind Ravichandran
#Sun Mar 12 00:41:05 CET 2017
#simulation.py
#Description:
 
from scipy import *
import numpy as np

class Simulation:
    """ data structure for storing general simulation information"""

    def __init__(self, lx, ly, dt, nsteps, nbeads, nsamp, nfils, nbpf, density, kappa, km, pa, bl, sigma):

        self.lx = lx
        self.ly = ly
        self.dt = dt
        self.nsteps = nsteps
        self.nbeads = nbeads
        self.nsamp = nsamp
        self.nfils = nfils
        self.nbpf = nbpf
        self.density = density
        self.kappa = kappa
        self.km = km
        self.pa = pa
        self.bl = bl
        self.sigma = sigma
        self.density = "{:.2f}".format(float(density))
        self.kappa = "{:.1f}".format(float(kappa))
        self.km = "{:.1f}".format(float(km))
        self.pa = "{:.1f}".format(float(pa))

        return
