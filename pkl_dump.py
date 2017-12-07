
""" Main script for Streaming (Levy Flight) type calculations """ 

import matplotlib as mpl
from read_write import read_data
#from matplotlib import pyplot as plt
import numpy as np
import pickle as pkl

from scipy.signal import correlate
from scipy.optimize import curve_fit

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
        
        data = fils.xi.transpose(2,1,0) 		# Adopt Guglielmo's convention here to avoid confusion
        pkl.dump(data, open('com.pkl','wb'))

        data = fils.ori.transpose(2,1,0) 		# Adopt Guglielmo's convention here to avoid confusion
        pkl.dump(data, open('ori.pkl','wb'))

    return
    
if __name__ == '__main__':
    main()    
