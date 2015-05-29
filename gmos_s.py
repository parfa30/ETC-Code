#!/usr/bin/env python
from numpy import interp
from scipy import interpolate
import scipy as sci

from ETC_inputs import *
from Instrument import *


class GMOS_S(Instrument):
    def __init__(self):
        Instrument.__init__(self,'/Users/parkerf/Documents/Research/ETC/ETC/inst_parfiles/gmos_s.txt')
        
        self.dI = self.par['dI']
        self.RON = self.par['RON']
        self.npix = self.par['pix_sampling'] #mapping of spaxel to pixels for IFUs. Usually 2 pixels/spax/wavelength. 1 for non IFUs
        self.pix = self.par['pix_size'] #pixel size in arcsec/pix
        self.slit_w = float(raw_input('GMOS_S slit width: (0.25, 0.5, 0.75, 1.0, 1.5, 2.0) '))
        self.grating = str(raw_input('Grating [R150, R400, R600, R831]: '))
        self.filter = str(raw_input('Filter [GG455,OG515,RG610,RG780]: '))
        R = {'R150':631,'R400':1918,'R600':3744,'R831':4396}
        disp = {'R150':0.193*10.,'R400':.074*10.,'R600':.052*10.,'R831':.038*10.} #A
        self.R = R[self.grating]
        self.disp = disp[self.grating]
        self.wave = np.arange(I.lower,I.upper, I.step)
        self.N_spec_pix = self.par['N_spec_pix'] #number pixels on detector along spectral direction
        self.center = I.center_input
        
    def spax(self):
        #size of spatial pixel in arcsec^2.
        return self.pix*self.pix

    def tel_tput(self):
        l = self.wave
        a=np.loadtxt('/Users/parkerf/Documents/Research/ETC/transmissions/aluminum_thr.txt')
        self.percent(a)
        degrade=a[:,1]*self.mdeg
        A=np.array(zip(a[:,0],degrade**self.nmirr))    
        Tel=interp(l,A[:,0],A[:,1])
        return Tel
        
    def n(self,R):
        #number of pixels per spectral pixel
        return np.rint(R*2./self.pix)

    def width_npix(self):
        return np.rint(self.slit_w/self.pix)

    def QE(self):
        data = np.genfromtxt('http://www.gemini.edu/sciops/instruments/gmos/gmos_n_ccd_hamamatsu_sc.txt')
        QE = np.interp(self.wave,data[:,0]*10,data[:,1])

        #QE(self.wave)
        return QE

    def Tput_slit(self):
        #This is a stand-in for gaussian aperture shape only
        
        tput = sci.special.erf((self.slit_w/2.)/(np.sqrt(2.)*I.FWHM/2.35482))
        return tput*np.ones(len(self.wave))

    def Tput_filter(self):
        #These are the spectroscopic blocking filters for GMOS_S.
        GG455 = '/Users/parkerf/Documents/Research/ETC/ETC/tput/gmos_s_gg455_G0329.txt'
        OG515 = '/Users/parkerf/Documents/Research/ETC/ETC/tput/gmos_s_og515_G0330.txt'
        RG610 = '/Users/parkerf/Documents/Research/ETC/ETC/tput/gmos_s_rg610_G0331.txt'
        RG780 = '/Users/parkerf/Documents/Research/ETC/ETC/tput/gmos_s_rg780_G0334.txt'
        
        filters = (GG455,OG515,RG610,RG780)
        filternames = ('GG455','OG515','RG610','RG780') 
        F = {}
        for i,f in enumerate(filters):
            data = np.genfromtxt(f)
            line = np.interp(self.wave,data[:,0]*10,data[:,1],left=0,right=0)
            F[filternames[i]] = line

        return F[self.filter]
    
    def Tput_grating(self):
        R150 = 'http://www.gemini.edu/sciops/instruments/gmos/gratings/gmos_n_R150_G5306.txt'
        R400 = 'http://www.gemini.edu/sciops/instruments/gmos/gratings/gmos_n_R400_G5305.txt'
        R600 = 'http://www.gemini.edu/sciops/instruments/gmos/gratings/gmos_n_R600_G5304.txt'
        R831 = 'http://www.gemini.edu/sciops/instruments/gmos/gratings/gmos_n_R831_G5302.txt'

        gratings = (R150,R400,R600,R831)
        gratingnames = ('R150','R400','R600','R831')
        G = {}
        for i,g in enumerate(gratings):
            data = np.genfromtxt(g)
            line = np.interp(self.wave,data[:,0]*10,data[:,1],left=0,right=0)
            G[gratingnames[i]] = line

        return G[self.grating]
        
    def Tput(self):
        #Throughput contains QE*filter*grating*slit. Telescope throughput is added in later. This is then used in sncosmo for flux of source. Need to multiply this by the BG though.
        
        Total = self.QE()*self.Tput_filter()*self.Tput_grating()#*self.Tput_slit()
        np.savetxt('/Users/parkerf/Documents/Research/ETC/ETC/gmos_s/'+self.filter+'.txt', np.transpose((self.wave,Total)),newline=os.linesep)
        return Total
      
    def plot_Tput(self):
        #I should rewrite this with a new program for just the outputs
        Tput = self.Tput()
        T = self.tel_tput()
        fig = py.figure(figsize=[7,5])
        ax  = fig.add_axes([0.15,0.1,0.8,0.8]) 
        ax.plot(self.wave,Tput*T) 
        py.title('Transmission of '+self.name+' with filter ' +self.filter+' and grating '+self.grating+' + Tel')
        py.ylabel('Throughput %')
        py.xlabel('Wavelength ($\AA$)')
        py.ylim(0,0.5)
        py.grid()
        py.minorticks_on()

    def plot_Tput_components(self):
        elements = (self.QE(),self.Tput_filter(),self.Tput_grating(),self.tel_tput(),self.Tput_slit())
        
        from itertools import cycle
        lines = ["-","--","-.",":","-"]
        linecycler = cycle(lines)
        fig = py.figure(figsize=[7,5])
        ax  = fig.add_axes([0.15,0.1,0.8,0.8])
        py.title('Elements of GMOS-S Transmission with '+self.filter+'filter and '+self.grating+'grating')
        py.ylabel('Throughput %')
        py.xlabel('Wavelength ($\AA$)')
        py.ylim(0,1.0)
        #py.xlim(3000,10000)
        for e in elements:
            py.plot(self.wave,e,next(linecycler))
      
        py.show() 


