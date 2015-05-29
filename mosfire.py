#!/usr/bin/env python
from numpy import interp
from scipy import interpolate

from ETC_inputs import *
from Instrument import *

class MOSFIRE(Instrument):
    def __init__(self):
        Instrument.__init__(self,'/Users/parkerf/Documents/Research/ETC/ETC/inst_parfiles/mosfire.txt')
        
        self.dI = self.par['dI']
        self.RON = self.par['RON']
        self.npix = self.par['pix_sampling'] #pixels per spaxel per spectral element
        self.wave = np.arange(I.lower,I.upper, I.step)
        self.slit_w = float(raw_input('slit width (arcsec): '))
        self.filter = str(raw_input('Band [all, H, J, K, Y]: '))
        R = {'H': 3660*(0.7/self.slit_w), 'J': 3310*(0.7/self.slit_w), 'K': 3620*(0.7/self.slit_w), 'Y': 3380*(0.7/self.slit_w)}
        disp = {'H': 1.629, 'J': 1.303, 'K': 2.17, 'Y': 1.086}
        center = {'H': 16500, 'J': 12500, 'K': 22000, 'Y': 10500}
        self.R = R[self.filter]
        self.disp = disp[self.filter]
        
        self.N_spec_pix = self.par['N_spec_pix']
        
        self.pix = self.par['pix_size'] #pixel size in arcsec/pix
        if I.Source == 'line':
            self.center = line
        else:
            self.center = center[self.filter]
            
    def spax(self):
        #size of pixel in arcsec^2
        return self.pix**2      

    def tel_tput(self):
        a=np.loadtxt('/Users/parkerf/Documents/Research/ETC/transmissions/aluminum_thr.txt')
        self.percent(a)
        degrade=a[:,1]*self.mdeg
        A=np.array(zip(a[:,0],degrade**self.nmirr))    
        Tel=interp(self.wave,A[:,0],A[:,1])
        return Tel
    
    def n(self,R):
        #number of pixels per spectral pixel
        return (R*2./self.pix)

    def width_npix(self):
        return (self.slit_w/self.pix)

    def Tput(self):
        
        # The throughput of the instrument includes the throughput of common optical elements, configurable elements, and the detector QE. Also includes telescope and atmosphere.

        h = 'http://www2.keck.hawaii.edu/inst/mosfire/data/throughput/measured/Heff_mosfire+tel.dat'
        j = 'http://www2.keck.hawaii.edu/inst/mosfire/data/throughput/measured/Jeff_mosfire+tel.dat'
        y = 'http://www2.keck.hawaii.edu/inst/mosfire/data/throughput/measured/Yeff_mosfire+tel.dat'
        k = 'http://www2.keck.hawaii.edu/inst/mosfire/data/throughput/measured/Keff_mosfire+tel.dat'
        
        bands = (h,j,y,k)
        Bandnames = ('H','J','Y','K')
        Tlist = {}

        for i,band in enumerate(bands): 
            data = np.genfromtxt(band)
            Tlist[Bandnames[i]] = np.interp(self.wave,data[:,0],data[:,1], left = 0, right = 0)
        
        Tlist['all'] = Tlist['H']+Tlist['J']+Tlist['Y']+Tlist['K'] 
        np.savetxt('/Users/parkerf/Documents/Research/ETC/ETC/mosfire/'+self.filter+'.txt', np.transpose((self.wave,Tlist[self.filter])),newline=os.linesep)
        return Tlist[self.filter]
      
    def plot_Tput(self):
        #I should rewrite this with a new program for just the outputs
        Tput = self.Tput()
        T = self.tel_tput()
        fig = py.figure(figsize=[7,5])
        ax  = fig.add_axes([0.15,0.1,0.8,0.8]) 
        ax.plot(self.wave,Tput*T) 
        py.title('Transmission of '+self.name+' + Tel')
        py.ylabel('Throughput %')
        py.xlabel('Wavelength ($\AA$)')
        py.ylim(0,1.0)
        py.grid()
        py.minorticks_on()


