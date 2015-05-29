# Using sncosmo

import sys, re, os
import numpy as np
import astropy
import sncosmo


# Pick cosmology
from astropy.cosmology import WMAP9 as cosmo
#from astropy.cosmology import FlatLambdaCDM
#cosmo = FlatLambdaCDM(H0=69.6, Om0=0.286)

from astropy import units as u
import astropy.constants as const
import matplotlib.pyplot as plt

# Parameters
z = 1.2    # Redshift of SN
obsfilter = 'f110w'   # This is the name of the sample band we use
wave = np.arange(8000,20000,10)  # This is a sample wavelength array for
                                # which to get a spectrum

# Other defined things
absmag_V = -19.3           # We use this to normalize (needed?)
magsystem = 'ab'           # The magnitude system (both ab and vega should work)
modelphase = 0             # The phase of the SN
template = 'salt2'         # The name of the SN template (e.g. salt2 or hsiao)

# When using new bands, assuing these are in newband directory
# e.g. from WFC3IR from http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/
bandinfo = {}
for filterfile in os.listdir( 'newband'):
    words = re.split('\.',filterfile)
    filterdata = np.genfromtxt( 'newband/'+filterfile )
    # Remove regions with zero throughput
    iNonzero = np.where( filterdata[:,2]>0 )[0]
    band = sncosmo.Bandpass(filterdata[iNonzero,1],filterdata[iNonzero,2],name=words[0])
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
model.set(x0=templatescale*10**(-0.4*DM.value),z=z) #this only works for salt2 model. need amplitude for hsiao model
# Derive the observed magnitude and flux in chosen filter
obsmag = model.bandmag(obsfilter,magsystem,[modelphase])
bandflux = model.bandflux(obsfilter, modelphase ) # Flux in photons/s/cm^2
print 'We get mag %s (%s) and flux %s photons/s/cm^2'%(obsmag,magsystem,bandflux)

# Get fluxes (in ergs/s/cm^2/A)
spectrum = model.flux(modelphase,wave)

plt.figure(1)
plt.plot(wave,spectrum)
plt.xlabel('Observer frame wavelength (A)')
plt.ylabel('Flux (ergs/s/sm^2/A)')
plt.savefig('test.png')

plt.show()
