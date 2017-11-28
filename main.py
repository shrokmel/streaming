
""" Main script for Streaming (Levy Flight) type calculations """ 

import matplotlib as mpl
from read_write import read_data
from matplotlib import pyplot as plt
import numpy as np
import pickle as pkl

from scipy.signal import correlate
from scipy.optimize import curve_fit

# Speed autocorrelation
def sacf(fils):
    data = fils.xi.transpose(2,1,0) 		# Adopt Guglielmo's convention here to avoid confusion

    

def vacf(fils):
    data = fils.xi.transpose(2,1,0) 		# Adopt Guglielmo's convention here to avoid confusion

# Just for plotting trajectories of MT COM
def vis(fils):

    rodids = range(4)				# just plot some rods

    f, ax = plt.subplots(1,1)
    ax.set_aspect('equal')

    for rodid in rodids:
        traj = cells.xi[5000:,:,rodid]		# just take latter half of sim
        traj_shift = traj - traj.mean(axis=0)
        ax.plot(traj_shift[:,0], traj_shift[:,1])	# subtract traj mean for cleanliness

# following normalisation doesn't really work
def step_length_hist_normed(n,fils):
    data = fils.xi.transpose(2,1,0) 		# Adopt Guglielmo's convention here to avoid confusion

    # data = [particle ID, (x,y), time] 
    step = data[:,:,n:] - data[:,:,:-n]		# step size dt = n
    drift= step.mean(axis=0)
    step = step - drift
    step_size = (step**2).sum(axis=1)**.5

    # p, bc and be are probability density, bin count and bin edges respectively
    p, be = np.histogram(distribution, bins=100, normed=True)	
    bc = (be[:-1]+be[1:])/2.

    ax1.loglog(bc,p)

    # normalised version
    step_size_mean = ((step**2).sum(axis=1).mean())**.5		# normalisation
    distribution = np.abs(step_size.flatten())/step_size_mean

    # p, bc and be are probability density, bin count and bin edges respectively
    p, be = np.histogram(distribution, bins=100, normed=True)	
    bc = (be[:-1]+be[1:])/2.
    ax2.loglog(bc,p)

def step_length_hist(n, fils):
    data = fils.xi.transpose(2,1,0) 		# Adopt Guglielmo's convention here to avoid confusion

    # data = [particle ID, (x,y), time] 
    step = data[:,:,n:] - data[:,:,:-n]		# step size dt = n
    drift= step.mean(axis=0)
    step = step - drift
    step_size = (step**2).sum(axis=1)**.5

    distribution = np.abs(step_size.flatten())

    # p, bc and be are probability density, bin count and bin edges respectively
    p, be = np.histogram(distribution, bins=100, normed=True)	
    bc = (be[:-1]+be[1:])/2.

    powerlaw    = lambda x, amp, index: amp*(x**index)
    exponential = lambda x, amp, index: amp*(index**x)

    long_time = 30
    popt, pcov   = curve_fit(powerlaw,    bc[long_time:], p[long_time:])
    popt2, pcov2 = curve_fit(exponential, bc[20:], p[20:])

    distribution = np.abs(step_size.flatten())/step_size_mean

    # p, bc and be are probability density, bin count and bin edges respectively
    p, be = np.histogram(distribution, bins=100, normed=True)	
    bc = (be[:-1]+be[1:])/2.

    step_size = (step2**2).sum(axis=1)**.5
    p, be = np.histogram(np.abs(step_size.flatten()), bins=100, normed=True)	
    bc = (be[:-1]+be[1:])/2.

    ax1.loglog(bc, p)
    ax1.loglog(bc[20:], powerlaw(bc[20:], *popt), c='k')

    ax2.semilogy(bc, p)
    ax2.semilogy(bc[10:], exponential(bc[10:], *popt2), c='k')

    ax1.set_ylabel('prob. density')
    ax1.set_xlabel('step length')

    ax2.set_ylabel('prob. density')
    ax2.set_xlabel('step length')
    #plt.xlim([0.08,1.5])

# This is a waste of time as Guglielmo rightly pointed out
def two_time_corr(ax,n, fils):
    data = fils.xi.transpose(2,1,0) 		# Adopt Guglielmo's convention here to avoid confusion

    # data = [particle ID, (x,y), time] 
    step = data[:,:,n:] - data[:,:,:-n]		# step size dt = n
    drift= step.mean(axis=0)
    step = step - drift
    step_size = (step**2).sum(axis=1)**.5

    x = step_size[:10,:-1].flatten()		# correlation between two groups of time steps
    y = step_size[:10,1:].flatten()
    ax.scatter(x,y,s=1) 

f, (ax1) = plt.subplots(1,1)
def main():

    # Loading any data in the 
    # phase space from the appropriate folder
    fname="out_fil.h5"
    base ="/work/jiff26/jiff2611/PROJECTS/effective/Jobs/170731-Active_Diffusion/rotation_data/aspect_ratio_10.0/"
    density=[0.4];	pa=[0.0,1.0];	pp=[0.0];
 
    folders=[] 
    for di in density:
        for pai in pa:
            for ppi in pp:
                folders.append(base+'density_'+str(di)+'/pa_'+str(pai)+'/pp_'+str(ppi)+'/')

    axes = (ax)
    # all functions here, iterate parameter set by set
    for folder in folders:
        fils, sim = read_data(folder, fname)
        #vis(fils) 		# visualise traj
        for ax,n in zip(axes,[1,10,100,1000]):
            two_time_corr(ax,n,fils)	# step length histogram (PD)
    plt.show() 
    #plt.savefig('step_length_loglog_semilog2.pdf')
    return
    
if __name__ == '__main__':
    main()    
