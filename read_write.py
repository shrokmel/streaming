#!/bin/python
 
#Copyleft Arvind Ravichandran
#Sat Mar 11 17:00:16 CET 2017
#read_write.py
#Description:

import argparse, os, h5py
import numpy as np
from simulation import Simulation
from cells import Cells

def read_data(folder, fname):
    """ read simulation data through hdf5 file"""

    ### access the file

    fpath = folder + fname
    assert os.path.exists(fpath), "out_fil.h5 does NOT exist for " + fpath
    fl = h5py.File(fpath, 'r')

    ### read in the positions of filaments
    
    xi  = np.array(fl['/positions/xi'], dtype=np.float32)       # center of mass of MT
    ori = np.array(fl['/positions/ori'], dtype=np.float32)      # <orientation> of MT

    ### read in the box info

    lx = fl['/info/box/x'][...]
    ly = fl['/info/box/y'][...]


    ### read in the general simulation info
    
    dt = fl['/info/dt'][...]
    nsteps = fl['/info/nsteps'][...]
    nbeads = fl['/info/nbeads'][...]
    nsamp = fl['/info/nsamp'][...]

    ### read in the filament information
    
    nfils = fl['/info/nfils'][...]
    nbpf = fl['/info/nbpf'][...]

    ### read in the simulation parameters

    density = fl['/param/density'][...]
    kappa = fl['/param/kappa'][...]
    km = fl['/param/km'][...]
    pa = fl['/param/pa'][...]
    bl = fl['/param/bl'][...]
    sigma = fl['/param/sigma'][...]

    ### close the file

    fl.close()

    sim = Simulation(lx, ly, dt, nsteps, nbeads, nsamp, nfils, nbpf, density, kappa, km, pa, bl, sigma)
    fils = Cells(xi, ori, sim)
    return fils, sim

def write_analysis_data(x, y, sim, fname):
    """ write analysis data to the corresponding file"""
    print(str(sim.density))

    base = "./../../DATA/"
    base += "DIS/lag5/"
    os.system("mkdir -p " + base)
    if fname:
        fpath = base + fname + ".txt"
    else:
        fpath = base + "spa1_density_" + \
            str(sim.density) + "_kappa_" + str(sim.kappa) + "_km_" + str(sim.km) + "_panti_" + str(sim.pa) + ".txt"

    fl = open(fpath, 'w')
    print(fpath)
    N = len(x)
    for j in range(N):
        fl.write(str(x[j]) + '\t\t' + str(y[j]) + '\n')

    fl.close()

    return
