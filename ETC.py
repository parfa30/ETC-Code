#!/usr/bin/env python
import os
import glob
import numpy as np
from scipy import integrate
from scipy import interpolate
from numpy import interp
import matplotlib.pyplot as py

from SNR import *

l = wave
cahk = np.where((l>np.rint((I.Redshift+1)*3925)) & (l<np.rint((I.Redshift+1)*3970)))
py.plot(l[cahk],Apt_SNR()[0][cahk])
py.show()  


'''
#Call the function 


from SNR import *






files = os.listdir('/Users/parkerf/Documents/Research/ETC/SN2011fe')
inputs = str(raw_input('single or multiple input files? (s/m)'))
if inputs == 'm':
    for file in glob.glob('/Users/parkerf/Documents/Research/ETC/SN2011fe/*.fit'):    
        t, z, FWHM, name, lambda_obs, flux, var = fits_file(file)
        wave = redshift(lambda_obs)
    
        #plot_Tput(inst,'dich_b')
        plot_SNR(inst)
        PSF_SNR()
        #SNR_from_file()
        #plot_BG(inst)
        
        py.savefig('/Users/parkerf/Documents/Research/ETC/'+str(inst.name)+'_ETC/snifs_'+str(name[1])+'.pdf')
        
elif inputs == 's':
    t, z, FWHM, name, lambda_obs, flux, var = fits_file('/Users/parkerf/Documents/Research/ETC/SN2011fe/11feP007.fit')
    wave = redshift(lambda_obs)
    
    
    #PSF_SNR()
    #plot_Tput(inst,'dich_b')
    #plot_BG(inst)
    #Ext_Atm(inst)
    
    #py.savefig('/Users/parkerf/Documents/Research/ETC/'+str(inst.name)+'_ETC/snifs_'+str(name[1])+'.pdf')
'''

