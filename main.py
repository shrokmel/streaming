
""" Main script for analysis calculations """ 

import matplotlib as mpl
from read_write import read_data
from matplotlib import pyplot as plt
import numpy as np
import pickle as pkl

from scipy.optimize import curve_fit

# Just for plotting
def vis(fname):
    base ="/work/jiff26/jiff2611/PROJECTS/effective/Jobs/170731-Active_Diffusion/rotation_data/aspect_ratio_10.0/"
    density=[0.4];	pa=[1.0];	pp=[0.0];
 
    folders=[] 
    for di in density:
        for pai in pa:
            for ppi in pp:
                folders.append(base+'density_'+str(di)+'/pa_'+str(pai)+'/pp_'+str(ppi)+'/')

    for folder in folders:
        cells, sim = read_data(folder, fname)

        rodids = range(4)

        f, ax1 = plt.subplots(1,1)

        for rodid in rodids:
            traj = cells.xi[5000:,:,rodid]

            #disp = traj - traj[0]

            #disp2 = np.sqrt(np.square(disp[:,0]) + np.square(disp[:,1]))

            ax1.plot(traj[:,0] - traj[:,0].mean(), traj[:,1] - traj[:,1].mean())
            ax1.set_aspect('equal')
            #ax2.plot(traj[:,0],traj[:,1])
            #ax2.set_aspect('equal')

def main():

    fname="out_fil.h5"
    vis(fname) 
    plt.show()
    
    return
    
if __name__ == '__main__':
    main()    
