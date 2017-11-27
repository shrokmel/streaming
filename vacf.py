
""" Main script for analysis calculations """ 

##############################################################################
import matplotlib as mpl
from read_write import read_data
from matplotlib import pyplot as plt
import numpy as np
import pickle as pkl

from scipy.optimize import curve_fit
       
##############################################################################

def vacf(fname):
    base ="/work/jiff26/jiff2611/PROJECTS/effective/Jobs/170731-Active_Diffusion/rotation_data/aspect_ratio_10.0/"
    density=[0.4];	pa=[0.0, 0.4, 0.8, 1.0];	pp=[0.0];
 
    folders=[] 
    for di in density:
        for pai in pa:
            for ppi in pp:
                folders.append(base+'density_'+str(di)+'/pa_'+str(pai)+'/pp_'+str(ppi)+'/')

    lag = 1

    plotter0 = np.zeros((50,2))
    plotter1 = np.zeros((50,2))

    for fi, folder in enumerate(folders):
        cells, sim = read_data(folder, fname)

        lags_v = [int(l) for l in np.logspace(0,np.log10(len(cells.xi)),50)]

        plotter = np.zeros((len(lags_v),2))

        traj = cells.xi[:,:,:]
        v = traj[:-lag] - traj[lag:]

        for i,lagv in enumerate(lags_v):
            v_corr = np.sum(v[:-lagv]*v[lagv:],axis=1)

            plotter[i,0] = lagv

            if len(v_corr)>0:
                plotter[i,1] = v_corr.mean()
            else:
                plotter[i,1] = 0 

        if fi==0:
            plotter0 = plotter 
            #plt.semilogx(plotter0[:,0], plotter0[:,1])
        else:
            plotter1 = plotter - plotter0
            plt.semilogy(lags_v, plotter1[:,1])

def disp(fname):
    base ="/work/jiff26/jiff2611/PROJECTS/effective/Jobs/170731-Active_Diffusion/rotation_data/aspect_ratio_10.0/"
    density=[0.4];	pa=[1.0];	pp=[0.0];
 
    folders=[] 
    for di in density:
        for pai in pa:
            for ppi in pp:
                folders.append(base+'density_'+str(di)+'/pa_'+str(pai)+'/pp_'+str(ppi)+'/')

    for folder in folders:

        cells, sim = read_data(folder, fname)
        rodids = [0,10]

        f, (ax1,ax2) = plt.subplots(1,2)

        for rodid in rodids:
            traj = cells.xi[5000:,:,rodid]

            print(traj.shape)

            disp = traj - traj[0]

            disp2 = np.sqrt(np.square(disp[:,0]) + np.square(disp[:,1]))

            ax1.plot(disp2)
            ax2.plot(traj[:,0],traj[:,1])
            ax2.set_aspect('equal')

def levy(fname):
    base ="/work/jiff26/jiff2611/PROJECTS/effective/Jobs/170731-Active_Diffusion/rotation_data/aspect_ratio_10.0/"
    density=[0.4];	pa=[0.2,1.0];	pp=[0.0];

    powerlaw = lambda x, amp, index: amp * (x**index)
 
    folders=[] 
    for di in density:
        for pai in pa:
            for ppi in pp:
                folders.append(base+'density_'+str(di)+'/pa_'+str(pai)+'/pp_'+str(ppi)+'/')

    for folder in folders:

        rodids = [153]
        #lags = (10, 100, 500)
        lags = (100,)
        long_time = 10
        for lag in lags:
            for rodid in rodids:
                cells, sim = read_data(folder, fname)
                traj = cells.xi[5000:,:,:]

                dr = traj[:-lag] - traj[lag:]
        
                disps = np.sqrt(np.square(dr[:,0]) + np.square(dr[:,1]))

            hist, be = np.histogram(disps/10,bins=100,normed=True)
            bw = be[1] - be[0]
            bc = be[1:] - bw/2

            popt, pcov = curve_fit(powerlaw, bc[long_time:], hist[long_time:])

            plt.loglog(bc, hist)
            plt.loglog(bc[long_time:], powerlaw(bc[long_time:],*popt),c='k')
            print(*popt)

def main():

    fname="out_fil.h5"
    vacf(fname) 
    plt.show()
    
    return
    
if __name__ == '__main__':
    main()    
