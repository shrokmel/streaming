
""" Main script for Streaming (Levy Flight) type calculations """ 

import matplotlib as mpl
from read_write import read_data
from matplotlib import pyplot as plt
import numpy as np
import pickle as pkl

from scipy.optimize import curve_fit

# Just for plotting trajectories of MT COM
def vis(fils):

    rodids = range(4)				# just plot some rods

    f, ax = plt.subplots(1,1)
    ax.set_aspect('equal')

    for rodid in rodids:
        traj = cells.xi[5000:,:,rodid]		# just take latter half of sim
        traj_shift = traj - traj.mean(axis=0)
        ax.plot(traj_shift[:,0], traj_shift[:,1])	# subtract traj mean for cleanliness

def step_length_hist(fils):
    data = fils.xi.transpose(2,1,0) 		# Adopt Guglielmo's convention here to avoid confusion

    # data = [particle ID, (x,y), time] 

    step = np.diff(data, axis=2)		# step size dt = 1 (ARBITRARY)

    # p, bc and be are probability density, bin count and bin edges respectively
    p, be = np.histogram(np.abs(step.flatten()), bins=100, normed=True)	
    bc = (be[:-1]+be[1:])/2.

    powerlaw = lambda x, amp, index: amp*(x**index)

    long_time = 30
    popt, pcov = curve_fit(powerlaw, bc[long_time:], p[long_time:])

    plt.loglog(bc, p)
    plt.loglog(bc[20:], powerlaw(bc[20:], *popt), c='k')

    plt.ylabel('prob. density')
    plt.xlabel('step length')
    #plt.xlim([0.08,1.5])

def main():

    # Loading any data in the 
    # phase space from the appropriate folder
    fname="out_fil.h5"
    base ="/work/jiff26/jiff2611/PROJECTS/effective/Jobs/170731-Active_Diffusion/rotation_data/aspect_ratio_10.0/"
    density=[0.4];	pa=[1.0];	pp=[0.0];
 
    folders=[] 
    for di in density:
        for pai in pa:
            for ppi in pp:
                folders.append(base+'density_'+str(di)+'/pa_'+str(pai)+'/pp_'+str(ppi)+'/')

    # all functions here, iterate parameter set by set
    for folder in folders:
        fils, sim = read_data(folder, fname)
        #vis(fils) 		# visualise traj
        step_length_hist(fils)	# step length histogram (PD)
   
    #plt.show() 
    plt.savefig('power_law.pdf')
    return
    
if __name__ == '__main__':
    main()    
