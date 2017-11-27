
""" Main script for analysis calculations """ 

##############################################################################
import matplotlib as mpl
from read_write import read_data
from matplotlib import pyplot as plt
import numpy as np
import pickle as pkl

from scipy.optimize import curve_fit
       
##############################################################################
def gug(fname):
    base ="/work/jiff26/jiff2611/PROJECTS/effective/Jobs/170731-Active_Diffusion/rotation_data/aspect_ratio_10.0/"
    density=[0.4];	pa=[1.0];	pp=[0.0];
 
    folders=[] 
    for di in density:
        for pai in pa:
            for ppi in pp:
                folders.append(base+'density_'+str(di)+'/pa_'+str(pai)+'/pp_'+str(ppi)+'/')

    for folder in folders:
        cells, sim = read_data(folder, fname)
        traj = cells.xi[:,:,:50]
        fp = open('./gug.pkl','wb')
        pkl.dump(traj,fp)

        fp.close()

def disp(fname):
    base ="/work/jiff26/jiff2611/PROJECTS/effective/Jobs/170731-Active_Diffusion/rotation_data/aspect_ratio_10.0/"
    density=[0.4];	pa=[1.0];	pp=[0.0];
 
    folders=[] 
    for di in density:
        for pai in pa:
            for ppi in pp:
                folders.append(base+'density_'+str(di)+'/pa_'+str(pai)+'/pp_'+str(ppi)+'/')

    for folder in folders:

        rodids = [0]

        f, (ax1,ax2,ax3) = plt.subplots(3,1)

        for rodid in rodids:
            cells, sim = read_data(folder, fname)
            traj = cells.xi[5000:,:,rodid]

            print(traj.shape)

            disp = traj - traj[0]

            disp2 = np.sqrt(np.square(disp[:,0]) + np.square(disp[:,1]))

            ax3.plot(traj[:1000,0],traj[:1000,1],c='b')
            ax3.plot(traj[1000:2000,0],traj[1000:2000,1],c='y')
            ax3.plot(traj[2000:3000,0],traj[2000:3000,1],c='g')
            ax3.plot(traj[3000:4000,0],traj[3000:4000,1],c='r')
            ax3.plot(traj[4000:5000,0],traj[4000:5000,1],c='k')

            ax3.set_aspect('equal')
       
        lag = 100;	thresh = 10

        slope = disp2[:-lag] - disp2[lag:]
        ax1.plot(np.abs(slope))

        #ax2.plot(np.abs(slope))
        ax2.plot(np.abs(slope)>thresh)

        plotter = slope

        plotter[np.abs(slope)>thresh] = 1
        plotter[np.abs(slope)<thresh] = 0
        #ax2.plot(plotter,'o')

def main():

    fname="out_fil.h5"
    gug(fname) 
    #plt.show()
    
    return
    
if __name__ == '__main__':
    main()    
