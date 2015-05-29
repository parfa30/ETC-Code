#!/usr/bin/env python
import pyfits
import sys, re, os
import numpy as np
import matplotlib.pyplot as py
import astropy
import sncosmo

############
###INPUTS###
############
#Inputs from a parameter file named 'ETC_inputs'#

class Inputs(object):
    def __init__(self,paramfile):

        pars = np.genfromtxt(paramfile, dtype = str)
        #print pars
        
        self.par = {}
        for p in pars:
            key = p[0]
            #print key
            value = p[1]
            try:
                value = float(value)
            except Exception as e:
                #print e
                pass
            self.par[key] = value
        
        self.Inst = self.par['Inst'] #mosfire, SNIFS, WFIRST
        self.t = self.par['t']
        self.N = self.par['N'] #number of exposures
        self.FWHM = self.par['FWHM']
        self.apt = self.par['aperture_shape'] #gaussian or moffat
        self.apt_r = self.par['aperture_rad'] #how large of an area you want for aperture photometry. Will be multiplied by FWHM. nominal is 3/4.
        self.amass = self.par['airmass'] #airmass
        self.reads = self.par['reads'] #Only using with MOSFIRE right now. probably can forget about for most part
        self.n_bin_spat = self.par['n_bin_spatial']
        
        self.Source = self.par['Source'] # salt2, hsiao, 11fe, twins, line
        self.Redshift = self.par['Redshift']
        
        self.lower = self.par['lower'] #lower end of wavelength range in Angstroms
        self.upper = self.par['upper'] #upper end of wavelength range in Angstroms
        self.step = self.par['step'] #delta_lambda. Not sure what best to put here in general. Is later convolved with dispersion.
        
        self.line = self.par['line'] #location of line flux in angstroms. can also be a wavelength that want the S/N for
        self.line_flux = self.par['line_flux'] #in ergs/s/cm^2/A
        self.FWHM_line = self.par['FWHM_line'] #in km/s

        self.option = self.par['option'] #Plots [a(SNR), b(SNR+Tput), c(SNR+Tput+BG)]
        self.windows = self.par['windows'] #sliced plots y/n
        self.save = self.par['save'] #do you want to save the plots or just show them: yes/no
        self.plot_inp = self.par['plot_input'] #do you want to show the input spectrum: yes/no

        self.loc = self.par['location'] #North:mk, South:cp
        self.wvap = self.par['wavervapor'] #MK: 1.0, 1.6, 3.0, 5.0mm; CP: 2.3, 4.3, 7.6, 10.0mm
        self.sky = self.par['sky'] #darkest,dark,grey,bright
        self.center_input = self.par['center'] #center wavelength. Sometimes this will be line, sometimes something else.

        
I = Inputs('/Users/parkerf/Desktop/ETC_inputs.txt')
z = I.Redshift
def redshift(lambda_obs):
    lambda_new = lambda_obs*(1.+ I.Redshift)
    return lambda_new
line = redshift(I.line)

#####################
###  Instrument   ###
#####################

from Instrument import *
from SNIFS import *
from mosfire import *
from gmos_s import *

if I.Inst.lower() == 'snifs':
    inst = SNIFS()
elif I.Inst.lower() == 'mosfire':
    inst = MOSFIRE()
elif I.Inst.lower() == 'gmos_s':
    inst = GMOS_S()
else:
    print "I don't have that instrument yet"

inst.Tput()
i_wave = np.arange(I.lower,I.upper,I.step)

#####################
###   SNCOSMO    ###
#####################

## This code is taken from SN cosmo/Jakob Nordin ##

def SNcosmo_template():
    # Pick cosmology
    from astropy.cosmology import WMAP9 as cosmo
    #from astropy.cosmology import FlatLambdaCDM
    #cosmo = FlatLambdaCDM(H0=69.6, Om0=0.286)

    from astropy import units as u
    import astropy.constants as const

    # Parameters
    obsfilter = inst.filter   # This is the name of the sample band we use

    # Other defined things
    absmag_V = -19.3           # We use this to normalize (needed?)
    magsystem = 'ab'           # The magnitude system (both ab and vega should work)
    modelphase = 0             # The phase of the SN
    template = I.Source        # The name of the SN template (e.g. salt2 or hsiao)

    # When using new bands, assuming these are in newband directory
    # e.g. from WFC3IR from http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/
    bandinfo = {}
    for filterfile in os.listdir(I.Inst):
        if filterfile == '.DS_Store':  
                pass
        else:
            words = re.split('\.',filterfile)
            filterdata = np.genfromtxt(I.Inst+'/'+filterfile )
            # Remove regions with zero throughput
            #iNonzero = np.where( filterdata[:,1]>0 )[0]
            band = sncosmo.Bandpass(i_wave,np.ones(len(i_wave)),name=words[0])
            #band = sncosmo.Bandpass(filterdata[:,0],filterdata[:,1],name=words[0])
            sncosmo.registry.register(band)

    # Scale template to chosen absolute V (Vega)
    model = sncosmo.Model(source=template)
    model.set(z=0)
    magdiff = model.bandmag('bessellv','vega',[0])-absmag_V
    templatescale = 10**(0.4*magdiff)
    print 'Get scale %s to match absolute mag.'%(templatescale)
    
    # Get distance modulus
    DM = cosmo.distmod(z)
    print 'Get Distance Modulus %s'%(DM)

    # Create a model assuming both of these scales
    model = sncosmo.Model(source=template)
    print model
    model.set(x0=templatescale*10**(-0.4*DM.value),z=z) #this only works for salt2 model. need amplitude for hsiao model
    
    # Derive the observed magnitude and flux in chosen filter
    obsmag = model.bandmag(obsfilter,magsystem,[modelphase])
    bandflux = model.bandflux(obsfilter, modelphase ) # Flux in photons/s/cm^2
    print 'We get mag %s (%s) and flux %s photons/s/cm^2'%(obsmag,magsystem,bandflux)

    return model.flux(modelphase,i_wave) #wave is in observer frame


#######################
###   FITS files    ###
#######################

#Unpack fits file and print out all info that comes with the flux. The input file might need to be updated if want to use these inputs correctly. This will likely only be done for SNIFS.
def fits_file(file):
        hdu = pyfits.open(file)
        #print hdu.info()
        #print hdu[0].data
        #print hdu[0].header
        exp_time = hdu[0].header['EXPTIME']
        airmass = hdu[0].header['AIRMASS']
        seeing = hdu[0].header['SEEING']
        name = hdu[0].header['OBJECT']#, str(hdu[0].header['TMAX'])
        
        crval = hdu[ 0 ].header[ "CRVAL1" ]
        cdelt = hdu[ 0 ].header[ "CDELT1" ]
        naxis = hdu[ 0 ].header[ "NAXIS1" ]
        lambda_obs = crval + cdelt * np.arange(naxis)
    
        flux = (hdu[0].data)
        var = (hdu[1].data)
        hdu.close()

        i_wave = redshift(lambda_obs)

        print "Make sure that the following are the same in your input file:"
        print 'exp_time: ' + str(exp_time)
        print 'airmass: ' + str(airmass)
        print 'seeing: ' + str(seeing)
        print 'name: ' + str(name)
        print 'wave: ' + str(i_wave)
        print 'var: ' + str(var)
        
        return flux

#######################
###   Line Flux    ###
#######################

## Test code for line fluxes to compare to MOSFIRE ETC ##
def Line(wave):
    ind = (np.abs(wave-line)).argmin()
    return ind




#######################
###   Get Fluxes    ###
#######################


def flux():
    if I.Source.lower() == 'salt2':  #  (in phtons/s/cm^2)
        flux_t1 = SNcosmo_template()
        flux_t2 = 'false'
    elif I.Source.lower() == 'salt2-extended': #  (in phtons/s/cm^2) 
        flux_t1 = SNcosmo_template()
        flux_t2 = 'false'
    elif I.Source.lower() == 'hsiao': #  (in phtons/s/cm^2)
        flux_t1 = SNcosmo_template()
        flux_t2 = 'false'
    elif I.Source.lower() == '11fe': #  (in ergs/s/cm^2/A)
        SNF_file = '/Users/parkerf/Documents/Research/ETC/SN2011fe/11feP007.fit'
        flux_t1 = 'false'
        flux_t2 = fits_file(SNF_file)
    elif I.Source.lower() == 'twin': #  (in ergs/s/cm^2/A)
        SNF_file = '/Users/parkerf/Documents/Research/ETC/SNF20080610-000/SNF20080610-000_P005765_merged.fits'
        flux_t1 = 'false'
        flux_t2 = fits_file(SNF_file)
    elif I.Source.lower() == 'line': #  (in ergs/s/cm^2/A)
        flux_t1 = 'false'
        flux_t2 = np.zeros(len(i_wave))
    else:
        print 'That input is not valid'

    #np.savetxt('/Users/parkerf/Documents/Research/ETC/ETC/inputspectrum.txt', np.transpose((wave/10.,flux)),newline=os.linesep)
    return flux_t1, flux_t2

def print_input():   
    if I.plot_inp == 'yes':
        py.plot(i_wave,flux)
        py.title('Input SN spectrum')
        py.ylabel('ergs/s/cm^2/A')
        py.xlabel('Angstroms')
        py.show()
    else:
        pass


print_input()



   
'''
Other stuff that I don't need right now but I might be able to use in the future.
    
np.savetxt('/Users/parkerf/Desktop/sn_file.txt',np.c_[wave*10.**-1,spectrum])
py.plot(wave*10**-1,spectrum)
py.show()    
def plot_input():
    fig = py.figure(figsize=[7,5])
    ax  = fig.add_axes([0.15,0.1,0.8,0.8]) 
    ax.plot(wave, spectrum )
    py.ylabel('Flux (ergs/s/sm^2/A)')
    py.xlabel('Observer frame wavelength (A)')
    py.title('Source flux from sncosmo')
    #py.xlim((Redshift+1)*3600,(Redshift+1)*4100) #CaH&K
    #py.savefig('test.png')
    
    py.show()

#plot_input()
'''

    






