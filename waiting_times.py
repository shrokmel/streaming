
""" Main script for Streaming (Levy Flight) type calculations """ 

import matplotlib as mpl
from read_write import read_data
from matplotlib import pyplot as plt
import numpy as np
import pickle as pkl

from scipy.ndimage.filters import median_filter
from scipy.ndimage.filters import gaussian_filter1d
from scipy.signal import correlate
from scipy.optimize import curve_fit

def lag_times_hist(sigma, threshold, p):
    running = []
    diffusing = []
    for j in range(50):
        mask = (plt.Normalize()(gaussian_filter1d(p[j],sigma=sigma))>threshold).astype(float)
        # and now we compute the lag times .. guess how?!
        true_to_false = np.diff(mask)
        where_is_it = np.where(true_to_false)
        lags = np.diff(where_is_it)[0]
        other = lambda x: {1:0,0:1}[x]
        assigned_lags = {mask[0]:lags[::2],other(mask[0]):lags[1::2]}
        running.extend(assigned_lags[1])
        diffusing.extend(assigned_lags[0])
    count,edges = np.histogram(running,bins=np.arange(0,5000,100));
    edges = (edges[1:]+edges[:-1])/2
    plt.semilogy(edges,count,label=str(sigma)+' '+str(threshold))
    plt.xlabel('lag times of running motion')
    plt.legend(loc='best')
    plt.grid(1)

def traj(fils,n=1,t=0.5,sigma=50):

    sigma = sigma
    threshold = t

    data = fils.xi.transpose(2,1,0) 		# Adopt Guglielmo's convention here to avoid confusion

    # data = [particle ID, (x,y), time] 
    step = data[:,:,n:] - data[:,:,:-n]		# step size dt = n
    step_size = (step**2).sum(axis=1)**.5

    # coloured wrt velocity
    #for j in range(3):
    #    ax1.scatter(data[j][0],data[j][1],s=.5,c=plt.cm.jet(plt.Normalize()(gaussian_filter1d(step_size[j],sigma=sigma))))

    d = step/step_size[:,np.newaxis,:]		# unit vector of velocity
    p = (d[:,:,n:]*d[:,:,:-n]).sum(axis=1)	# vel correlation between n steps

    # coloured wrt orientation corr.
    #for j in range(3):
    #    ax2.scatter(data[j][0],data[j][1],s=.5,c=plt.cm.jet(plt.Normalize()(gaussian_filter1d(p[j],sigma=sigma))))

    for j in range(3):
        ax2.scatter(data[j][0],data[j][1],s=.5,c=plt.cm.jet(plt.Normalize()(gaussian_filter1d(step_size[j],sigma=sigma))))
        #ax2.scatter(data[j][0],data[j][1],s=.5)

        #threshold = 0.3
        mask = (plt.Normalize()(gaussian_filter1d(step_size[j],sigma=sigma))>threshold).astype(float)
        ax1.scatter(data[j][0],data[j][1],s=.5,c=plt.cm.gray_r(mask))

        #threshold = 0.5
        #mask = (plt.Normalize()(gaussian_filter1d(p[j],sigma=sigma))>threshold).astype(float)
        #ax2.scatter(data[j][0],data[j][1],s=.5,c=plt.cm.gray_r(mask))

    return

def traj_distance(fils,n=1,t=0.4,sigma=50):
    data = fils.xi.transpose(2,1,0) 		# Adopt Guglielmo's convention here to avoid confusion

    # data = [particle ID, (x,y), time] 
    step = data[:,:,n:] - data[:,:,:-n]		# step size dt = n
    step_size = (step**2).sum(axis=1)**.5

    #d = step/step_size[:,np.newaxis,:]		# unit vector of velocity
    #p = (d[:,:,n:]*d[:,:,:-n]).sum(axis=1)	# vel correlation between n steps

    #j = 1

    threshold = t
    #mask1 = (plt.Normalize()(gaussian_filter1d(step_size[j],sigma=sigma))>threshold).astype(float)[1:]
    #ax1.plot(mask1)

    #threshold = 0.6
    #mask2 = (plt.Normalize()(gaussian_filter1d(p[j],sigma=sigma))>threshold).astype(float)
    #ax2.plot(mask2)

    running = []
    diffusing = []

    for j in range(1250):     
        mask = (plt.Normalize()(gaussian_filter1d(step_size[j],sigma=sigma))>threshold).astype(float)[1:]
        true_to_false = np.diff(mask)

        # DISTANCE covered
        d = data[j,:,1:-2]		# prune to accomodate true_to_false dimensions
        stream_enter = d[:,true_to_false== 1]
        stream_exit  = d[:,true_to_false==-1]

        # if filament starts diffusive & ends diffusive
        if mask[0] == 0 and mask[-1] == 0:
            # streaming
            dist = stream_exit - stream_enter
            dist = (dist**2).sum(axis=0)**.5

            for ll in dist:
                running.append(ll)

            # diffusing
            dist = stream_enter[:,1:] - stream_exit[:,:-1]
            dist = (dist**2).sum(axis=0)**.5

            for ll in dist:
                diffusing.append(ll)

        # if filament starts streaming & ends diffusive
        elif mask[0] == 1 and mask[-1] == 0:
            # streaming
            dist = stream_exit[:,1:] - stream_enter
            dist = (dist**2).sum(axis=0)**.5
            for ll in dist:
                running.append(ll)

            # diffusing
            dist = stream_enter - stream_exit[:,:-1]
            dist = (dist**2).sum(axis=0)**.5

            for ll in dist:
                diffusing.append(ll)

        # if filament starts streaming & ends streaming
        elif mask[0] == 1 and mask[-1] == 1:
            # streaming
            dist = stream_exit[:,1:] - stream_enter[:,:-1]
            dist = (dist**2).sum(axis=0)**.5
            for ll in dist:
                running.append(ll)
            # diffusing
            dist = stream_enter - stream_exit
            dist = (dist**2).sum(axis=0)**.5
            for ll in dist:
                diffusing.append(ll)
        
        # if filament starts diffusive & ends streaming
        elif mask[0] == 0 and mask[-1] == 1:
            # streaming
            lags = stream_exit - stream_enter[:,:-1]
            dist = (dist**2).sum(axis=0)**.5
 
            if type(dist) is np.float64:
                running.append(dist)
            else:
                print(type(dist))
                for ll in dist:
                    running.append(ll)

            # diffusing
            dist = stream_enter[:,1:] - stream_exit
            dist = (dist**2).sum(axis=0)**.5
            for ll in dist:
                diffusing.append(ll)
    
    p, be = np.histogram(running, bins=100)
    bc = (be[:-1]+be[1:])/2.
    ax3.semilogy(bc,p,'o',label=str('thresh='+str(t)))
    ax3.set_xlabel('distance covered streaming')


    p, be = np.histogram(diffusing, bins=100)
    bc = (be[:-1]+be[1:])/2.

    ax4.semilogy(bc,p,'o',label=str('thresh='+str(t)))
    ax4.set_xlabel('distance covered in bundle')
    plt.legend()

    
    return

def traj_time(fils,n=1,t=0.1,sigma=30):
    data = fils.xi.transpose(2,1,0) 		# Adopt Guglielmo's convention here to avoid confusion

    # data = [particle ID, (x,y), time] 
    step = data[:,:,n:] - data[:,:,:-n]		# step size dt = n
    step_size = (step**2).sum(axis=1)**.5

    d = step/step_size[:,np.newaxis,:]		# unit vector of velocity
    p = (d[:,:,n:]*d[:,:,:-n]).sum(axis=1)	# vel correlation between n steps

    #j = 1

    threshold = t
    #mask1 = (plt.Normalize()(gaussian_filter1d(step_size[j],sigma=sigma))>threshold).astype(float)[1:]
    #ax1.plot(mask1)

    #threshold = 0.6
    #mask2 = (plt.Normalize()(gaussian_filter1d(p[j],sigma=sigma))>threshold).astype(float)
    #ax2.plot(mask2)

    running = []
    diffusing = []

    for j in range(1250):     
        mask = (plt.Normalize()(gaussian_filter1d(step_size[j],sigma=sigma))>threshold).astype(float)[1:]
        true_to_false = np.diff(mask)

        # TIME
        stream_enter = np.where(true_to_false == 1)[0]
        stream_exit  = np.where(true_to_false ==-1)[0]

        # if filament starts diffusive & ends diffusive
        if mask[0] == 0 and mask[-1] == 0:
            # streaming
            lags = stream_exit - stream_enter
            for ll in lags:
                running.append(ll)

            # diffusing
            lags = stream_enter[1:] - stream_exit[:-1]
            for ll in lags:
                diffusing.append(ll)

        # if filament starts streaming & ends diffusive
        elif mask[0] == 1 and mask[-1] == 0:
            # streaming
            lags = stream_exit[1:] - stream_enter
            for ll in lags:
                running.append(ll)

            # diffusing
            lags = stream_enter - stream_exit[:-1]
            for ll in lags:
                diffusing.append(ll)

        # if filament starts streaming & ends streaming
        elif mask[0] == 1 and mask[-1] == 1:
            # streaming
            lags = stream_exit[1:] - stream_enter[:-1]
            for ll in lags:
                running.append(ll)

            # diffusing
            lags = stream_enter - stream_exit
            for ll in lags:
                diffusing.append(ll)
        
        # if filament starts diffusive & ends streaming
        elif mask[0] == 0 and mask[-1] == 1:
            # streaming
            lags = stream_exit - stream_enter[:-1]
            for ll in lags:
                running.append(ll)

            # diffusing
            lags = stream_enter[1:] - stream_exit
            for ll in lags:
                diffusing.append(ll)
    
    p, be = np.histogram(running, bins=100)
    bc = (be[:-1]+be[1:])/2.
    ax5.semilogy(bc,p,'o',label=str('sigma='+str(sigma)))
    ax5.set_xlabel('time spent streaming')


    p, be = np.histogram(diffusing, bins=100)
    bc = (be[:-1]+be[1:])/2.

    ax6.semilogy(bc,p,'o',label=str('sigma='+str(sigma)))
    ax6.set_xlabel('time spent in bundle')
    plt.legend()

    #p, be = np.histogram(diffusing, bins=50)
    #bc = (be[:-1]+be[1:])/2.

    #plt.plot(bc,p)
    

    #plt.hist(diffusing,bins=50)    
    #plt.hist(running,bins=50)    
 
    '''
    # test if vel corr and step length correspond for threshold pair
   
    a = np.zeros((3,121))
    ct=0
    for ti in np.linspace(0,1,11):
        for tj in np.linspace(0,1,11):
            mask_step = (plt.Normalize()(gaussian_filter1d(step_size[j],sigma=sigma))>ti).astype(float)[1:]
            mask_velc = (plt.Normalize()(gaussian_filter1d(p[j],sigma=sigma))>tj).astype(float)
            a[0,ct] = ti
            a[1,ct] = tj
            a[2,ct] = np.sum(mask_step*mask_velc)
            ct = ct+1

    print(np.max(a[2,:]))
    print(np.min(a[2,:]))

    plt.scatter(a[0,:], a[1,:], c=a[2,:], cmap="cool")
    '''

    return

f, (ax1,ax2) = plt.subplots(1,2, figsize=(20,10))
ax1.set_aspect('equal')
ax2.set_aspect('equal')

f2, ((ax3,ax4),(ax5,ax6)) = plt.subplots(2,2) #figsize=(20,10))

ax3.set_xlim([0,200])
ax4.set_xlim([0,200])

ax5.set_xlim([0,8000])
ax6.set_xlim([0,8000 ])
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
    
        for t in [0.3]:
             #traj_distance(fils,t=t)        
             #traj_time(fils,t=t)
             traj(fils,t=t,sigma=30)
             traj_distance(fils,t=t, sigma=30)	
             traj_time(fils,t=t, sigma=30)        
    plt.show()
    return
    
if __name__ == '__main__':
    main()    
