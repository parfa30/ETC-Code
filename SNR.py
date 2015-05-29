#!/usr/bin/env python
import glob
from scipy import integrate
from scipy import interpolate
from numpy import interp

from ETC_inputs import *
from bg import *



###########################
#           SNR          #
########################## 
def source_flux():
    #Source is characterized in ETC_inputs. Flux comes as photons or energy/s/cm^2. Run through telescope and atmosphere. It goes through the instrument in sncosmo.
    T = inst.tel_tput()*Ext_Atm()*inst.apt()  #cm^2
    if flux()[1]== 'false':
        Flux = flux()[0]*T
    elif flux()[0]=='false':
        Flux = flux()[1]*T*((i_wave)/(h*c))*inst.Tput()
        
    return Flux #photons/sec/A
    
def res_convol(input_flux,input_wave):
    #Following MOSFIRE code, which does the convolution in velocity space

    vel = np.arange(1,200001,1)-100000
    wave2vel = (input_wave/inst.center-1.)*c*10**(-13) #cm in km/s
    interp_flux = np.interp(vel,wave2vel,input_flux)

    sigma = (c*10**(-13)/inst.R)*(1/(2*np.sqrt(2*np.log(2))))
    vel_kernel = np.arange(0,np.rint(8*sigma),1.)-np.rint(8*sigma)/2.
    gauss_kernel = (1/(np.sqrt(2*np.pi)*sigma))*np.exp(-0.5*(vel**2)/(sigma**2))

    convol_flux = np.convolve(interp_flux,gauss_kernel,'same')
    convol_wave = inst.center*((vel/(c*10**(-13)))+1)

    real_wave = np.arange(1.,inst.N_spec_pix+1.,1.)*inst.disp
    real_wave = real_wave - real_wave[np.rint(inst.N_spec_pix/2.)] + inst.center

    real_flux = np.interp(real_wave,convol_wave,convol_flux)

    return real_flux, real_wave
    
def line_flux(x_wave):
    #width of the line
    i_width = inst.center*I.FWHM_line*c**(-1)*10**(-13)
    true_width = np.sqrt(i_width**2+(inst.center/inst.R)**2.)
    sigma = true_width/(2*np.sqrt(2*np.log(2)))
    gauss_dist = 1/(np.sqrt(2*np.pi)*sigma)*np.exp(-0.5*(x_wave-inst.center)**2./(sigma**2))

    #flux of line
    line_ind = np.where(np.abs(x_wave-inst.center) <= 0.5*true_width)
    F = np.zeros(len(x_wave))
    F[line_ind] = I.line_flux
    
    return F*gauss_dist, line_ind, true_width

def line_SNR():
    t, Nexp, spax,pix, reads = I.t, np.sqrt(I.N), inst.spax(), inst.pix, I.reads

    #convolve throughput with mosfire resolution

    atm = res_convol(Ext_Atm(),i_wave)[0]
    tput = res_convol(inst.Tput(),i_wave)[0]#inst.tel_tput()*

    #convolve background and call the line flux
    Wave = res_convol(source_flux(),i_wave)[1]
    source = line_flux(Wave)[0]*((Wave)/(h*c))*inst.apt()*inst.disp*tput #*(I.line/inst.R)*atm
    bckgd = res_convol(BG()[0],i_wave)[0]
    line_ind = line_flux(Wave)[1]

    #area of enclosed energy
    R_in = (I.FWHM)/2. #apt.r should usually be 1.4 - check this is the case. 
    nspec = line_flux(Wave)[2]/inst.disp #npix in spectral direction
    nspat = inst.n(R_in) # number of pixels in spatial direction

    signal = source*t
    dI = inst.dI*np.ones(len(Wave)) #Dark Current
    RON = inst.RON*np.ones(len(Wave)) #Read Out Noise
    bckgd = bckgd*t*nspat*nspec
    noise = np.sqrt(signal+bckgd+nspat*dI*t+((nspat*RON**2)/I.reads))

    SNR = Nexp*(signal/noise)

    print "### print statements ###"
    print "Tput = " + str(tput[line_ind])
    print "line wavelength = " +str(np.median(Wave[line_ind]))
    print "Resolution = " +str(inst.center/inst.R)
    print "length along aperture = " + str(R_in*2.)
    print "number pixels in aperture = " +str(nspat)
    print "number pixels in dispersion direction = " +str(nspec)
    print "dark current per FWHM = " +str(np.median(nspat*nspec*(dI*t)[line_ind]))
    print "RON per FWHM = " +str(np.median(np.sqrt((nspec*nspat)/I.reads)*RON[line_ind]))
    print "background per FWHM = " +str(np.sum(bckgd[line_ind])*nspec)#str(np.mean(bckgd[line_ind])*nspec) 
    print "Signal per FWHM = " + str(np.sum(signal[line_ind]))#str(np.sum(signal[line_ind]))#
    print "Noise per FWHM = " + str(np.sum(np.sqrt(signal[line_ind]+bckgd[line_ind]*nspec+nspec*nspat*(RON**2)[line_ind]/I.reads+nspec*nspat*(dI*t)[line_ind])))#str(np.mean(noise[line_ind])*np.sqrt(nspec))
    print 'SNR per FWHM = ' + str(np.sum(signal[line_ind]/np.sqrt(signal[line_ind]+bckgd[line_ind]*nspec+nspec*nspat*(RON**2)[line_ind]/I.reads+nspec*nspat*(dI*t)[line_ind])))#str(np.mean(SNR[line_ind])*np.sqrt(nspec))#str(np.mean(SNR[line_ind]*np.sq]t(nspec)))
    print "##################"
    

    return SNR, signal, noise, Wave




def f(x,y):
    alpha = I.FWHM/2.35
    return ((alpha**2)*(2*np.pi))**(-1.)*np.exp(-(1./2.)*(x/alpha)**2.)*np.exp(-(1./2.)*(y/alpha)**2.)

def Gaussian(xlim,ylim):
    #G_area = integrate.quad(lambda r: ((alpha**2)*(2*np.pi))**(-1.)*np.exp(-(1./2.)*(r/alpha)**2.)*2*np.pi*r,0,limit)
    G_area = integrate.nquad(f, [[0,xlim],[0,ylim]])
    #G_area = integrate.quad(lambda r: ((beta-1)/(np.pi*Alpha**2))*(1+(r/Alpha)**2)**(-beta)*2.*np.pi*r,0,limit)
    return G_area

def Moffat(limit):
    beta = 3.
    alpha = I.FWHM/(2.*np.sqrt(2**(1./beta)-1.))
    M_area = integrate.quad(lambda r: ((alpha**2)*(2*np.pi))**(-1.)*np.exp(-(1./2.)*(r/alpha)**2.)*2*np.pi*r,0,limit)
    return M_area


    
def Apt_SNR():
    t, Nexp, spax,pix, reads = I.t, np.sqrt(I.N), inst.spax(), inst.pix, I.reads
    

    #area of enclosed energy
    R_in = (I.FWHM*I.apt_r)/2.
    nw = inst.width_npix()
    ylim = (nw*inst.pix)/2.
    n = inst.n(R_in) # number of pixels in spatial direction
    R_round = n * inst.pix/2.
    #area = n*pix
    n_bin_spat = n/I.n_bin_spat

    #convolve fluxes with resolution
    source = res_convol(source_flux(),i_wave)[0]
    bckgd = res_convol(BG()[0],i_wave)[0]
    wave = res_convol(source_flux(),i_wave)[1]
 
    line_ind = Line(wave)
     
    #select aperture type
    if I.apt == 'gaussian':
        apt = Gaussian(R_round,ylim)[0]
        tot = Gaussian(np.inf,np.inf)[0]
    elif I.apt == 'moffat':
        apt = Moffat(R_round)[0]
        tot = Moffat(np.inf)[0]
    else:
        print 'Please select an aperture shape (gaussian or moffat)'

    enclosed_energy = apt/tot

    signal = source*t*enclosed_energy
    
    dI = inst.dI*np.ones(len(wave)) #Dark Current
    RON = inst.RON*np.ones(len(wave)) #Read Out Noise
    bckgd = bckgd*t*n*nw
    noise = np.sqrt(signal+bckgd+n*dI*t+n*RON**2)

    SNR = Nexp*(signal/noise)

    print "### print statements ###"
    print "line wavelength = " +str(wave[line_ind])
    print "Resolution = " +str(center/inst.R)
    print "energy enclosed = " + str(enclosed_energy)
    print "length along aperture = " + str(R_round*2.)
    print "number pixels in aperture = " +str(n)
    print "number pixels in slit width = " +str(nw)

    print "dark current = " +str(n*(dI*t)[line_ind])
    print "RON = " +str((n*RON)[line_ind])
    print "bckd counts = " +str(bckgd[line_ind]/n) #working on something
    print "sqrt bckgd = " +str(np.sqrt(bckgd)[line_ind])
    print "Signal = " + str(signal[line_ind])
    print "Noise = " + str(noise[line_ind])
    print 'sqrt of signal = ' + str(np.sqrt(signal[line_ind]))
    print 'SNR = ' +str(SNR[line_ind])
    #print 'area/apt: ' +str(area/apt)
    print "##################"

    return SNR, signal, noise, wave

def noise():
    #Not sure what to do with this anymore
    noise = Apt_SNR()[2]
    gauss = np.random.normal(0,1,len(noise))
    #count, bins, ignored = py.hist(gauss, 30, normed=True)
    #py.plot(bins, 1/(1 * np.sqrt(2 * np.pi)) * np.exp( - (bins)**2 / (2 * 1**2)), linewidth=2, color='r')
    noise_A = noise * gauss
    py.plot(Apt_SNR()[3],gauss)
    py.show()
    return noise_A

def noisy_spectrum():
    return Apt_SNR()[1]+Apt_SNR()[2]

width = 2.
def binned_SNR():
    #new_spectrum = noise_spectrum() + Aperture()[0]
    noisy_data = Aperture()[1]+noise()
    binned_data = noisy_data.reshape(-1, width).mean(axis=1)
    binned_noise = np.std(noisy_data)
    #binned_SNR = binned_data/binned_noise
    #binned_SNR = binned_data/Aperture()[2]

    #digitized = np.digitize(new_spectrum, bins)
    #bin_means = [new_spectrum[digitized == i].mean() for i in range(1, len(bins))]
    return noisy_data, binned_data



## def new_data():
##     new_wave = np.arange(7000.,19988.+6.,12.)
##     New_data = interp(wave,new_wave,binned_SNR()[1])
##     py.plot(wave,Aperture()[1],wave,New_data)
##     py.plot(new_wave,binned_SNR()[1])
##     np.sqrt(New_data))
##     py.xlim((Redshift+1)*3800,(Redshift+1)*4100)
##     f, (ax1,ax2) = plt.subplots(2)
##     py.title('Template signal vs noisy signal')
##     ax1.plot(new_wave,binned_SNR()[1])
##     ax1.set_xlim((Redshift+1)*3800,(Redshift+1)*4100)
##     ax2.plot(wave,Aperture()[1])
##     ax2.set_xlim((Redshift+1)*3800,(Redshift+1)*4100)

##     py.show()


##     f, (ax1, ax2) = plt.subplots(2)
##     ax1.plot(new_wave,binned_SNR()[1]/np.sqrt(binned_SNR()[1]))
##     ax1.set_xlim((Redshift+1)*3800,(Redshift+1)*4100)
##     ax2.plot(wave,Aperture()[0])
##     ax2.set_xlim((Redshift+1)*3800,(Redshift+1)*4100)
##     py.title('SNR comparison')
##     py.show()


#print 'noisy stddev ' + str(np.std(binned_SNR()[0], axis = 0))    
#print 'binned stddev ' + str(np.std(binned_SNR()[1], axis = 0))         
        

###########################
#       OUTPUTS          #
########################## 

def plot_full():
    Source, Redshift, t = I.Source, I.Redshift, I.t
    
    if Source == 'line':
        fig, (ax1) = py.subplots(1,figsize=[7,5])
        ax1.plot(line_SNR()[3],line_SNR()[0], 'b-')
        ax1.set_ylabel('SNR per spectral pixel')
        #ax2.plot(o_wave,noisy_spectrum(), 'r')
        #ax2.set_ylabel('e- per spectral pixel')
        ax1.set_title('SNR and Noisy Spectrum from '+Source+ ' with ' + inst.name + ' @ z='+str(Redshift) + ' w/ t='+str(t),fontsize = 10 )
    else:
        fig, (ax1, ax2) = py.subplots(2, sharex =True, figsize=[7,5])
        ax1.plot(Apt_SNR()[3],Apt_SNR()[0], 'b-')
        ax1.set_ylabel('SNR per spectral pixel')
        ax2.plot(Apt_SNR()[3],noisy_spectrum(), 'r')
        ax2.set_ylabel('e- per spectral pixel')
        ax1.set_title('SNR and Noisy Spectrum from '+Source+ ' with ' + inst.name + ' @ z='+str(Redshift) + ' w/ t='+str(t),fontsize = 10 )
    
def plot_windows():
    Source, Redshift, t = I.Source, I.Redshift, I.t
    fig, (ax1, ax2, ax3, ax4, ax5) = py.subplots(5, figsize = [7,5])
    l = o_wave
    #CaHK
    cahk = np.where((l>np.rint((Redshift+1)*3925)) & (l<np.rint((Redshift+1)*3970)))
    ax1.plot(l[cahk],noisy_spectrum()[cahk])
    ax1.set_ylabel('CaHK')
    #SiII
    si2 = np.where((l>np.rint((Redshift+1)*4120)) & (l<np.rint((Redshift+1)*4165)))
    ax2.plot(l[si2],noisy_spectrum()[si2])
    ax2.set_ylabel('SiII')
    #SII
    s2 = np.where((l>np.rint((Redshift+1)*5300)) & (l<np.rint((Redshift+1)*5600)))
    ax3.plot(l[s2],noisy_spectrum()[s2])
    ax3.set_ylabel('SII')
    #O2
    o2 = np.where((l>np.rint((Redshift+1)*5000)) & (l<np.rint((Redshift+1)*5020)))
    ax4.plot(l[o2],noisy_spectrum()[o2])
    ax4.set_ylabel('O2')
    #Ha
    Ha = np.where((l>np.rint((Redshift+1)*6555)) & (l<np.rint((Redshift+1)*6570)))
    ax5.plot(l[Ha],noisy_spectrum()[Ha])
    ax5.set_ylabel('Halpha')

    ax1.set_title('Noisy Spectrum '+Source+ ' of ' + inst.name + ' @ z='+str(Redshift) +  'w/ t='+str(t) + ' at SN features')

    if I.windows == 'yes':
        py.show()
    elif I.windows == 'no':
        pass
    else:
        pass

def call_function():
    #Plots only what you want to see. Seen inputs for options.
    if I.option == 'a':
        if I.save == 'no':
            plot_full()
            py.show()
        elif I.save == 'yes':
            py.savefig('snr')
        
    elif I.option == 'b':
        if I.save == 'no':
            plot_full()
            inst.plot_Tput()
            py.show()
        elif I.option == 'yes':
            py.savefig(pdf,'source_'+str(t))
            py.savefig(pdf,'tput')
        
    elif I.option == 'c':
        if I.save == 'no':
            plot_full()
            inst.plot_Tput()
            plot_BG()
            py.show()
        elif I.save == 'yes':
            py.savefig(pdf,'source_'+str(t))
            py.savefig(pdf,'tput')
            py.savefig(pdf,'BG')
    else:
        print 'pick an option'      
        
call_function()
#plot_windows()
#inst.plot_Tput_components()


## sys.exit()
## def plot_SNR():
##     l = inst.L

##     apt=Aperture()
##     psf = PSF_SNR()
##     from_file = SNR_from_file()
    
##     if inst.name == 'SNIFS':
##         if plots.lower() == 'a':
##             f, (ax1) = py.subplots(1, sharex = True)
##             ax1.set_title('SNR for '+str(name)+ ' with ' + inst.name )
##             #ax1.plot(l,from_file)
##             #ax1.set_ylabel('Measured')
##             ax1.plot(l,apt)
##             ax1.set_ylabel('Aperture Photometry')
##             py.xlim(l[0],l.max())
            
##         elif plots.lower() == 'b':
##             f, (ax1, ax2) = py.subplots(2, sharex = True)
##             ax1.plot(l,apt)
##             ax1.set_ylabel('Aperture Photometry')
##             ax2.plot(l, psf)
##             ax2.set_ylabel('PSF Photometry')
##             py.xlim(l[0],l.max())
##         else:
##             print 'what PSF'
            
            
##     else:
##         if plots.lower() == 'a':
##             py.figure()
##             py.title('SNR for '+str(name)+ ' with ' + inst.name )
##             py.plot(l,apt)
##             py.ylabel('Aperture Photometry')   
            
##         elif plots.lower() == 'b':
##             f, (ax1, ax2) = py.subplots(2, sharex = True)
##             ax1.set_title('SNR for '+str(name)+ ' with ' + inst.name  )
##             ax1.plot(l,apt)
##             ax1.set_ylabel('Aperture Photometry')
##             ax2.plot(l,psf)
##             ax2.set_ylabel('PSF Photometry')
##         else:
##             print 'what PSF?'


##         if np.amax(np.sqrt(signal)/np.sqrt(signal+bckgd+inst.n(I.FWHM)*(inst.dI*t+inst.RON**2/reads))) >= 0.5:
##         print "This observation is photon noise limited."
##     else:
##         print "This observation is background limited."

'''
def SNR_from_file(inst):
    #SNR calculated from signal and variance from fits file. noise=sqrt(var)
    #This should be compared to the SNR calculated below with the model for noise and throughput.
    
    signal = source_flux()[0]*I.t
    noise = source_flux()[1]*I.t
    SN = signal/noise
    
    return SN

def PSF_pix(x,y):
    #Pixelizing the PSF for PSF photometry#
    xx, x0, yy, y0 = inst.grid()
    
    beta = 3.
    alpha = I.FWHM/(2.*np.sqrt(2**(1./beta)-1.))
    Alpha = I.FWHM/2.
    
    G = (2.*np.pi*Alpha)**(-2.)*np.exp(-((x-x0)**2.+(y-y0)**2.)/(2.*(Alpha**2.))) #
    M = (2.*np.pi)**(-1.)*(beta-1.)*(np.pi*alpha**2.)**(-1.)*(1+((x-x0)**2.+(y-y0)**2.)/alpha**2.)**(-beta) #

    return G, M

def PSF_SNR():
    #PSF photometry#
    xx, x0, yy, y0 = inst.grid()
    x, y = np.meshgrid(xx, yy) 
    t = I.t
    
    SNR = [] 
    P = PSF_pix(x,y)[1]
    #py.pcolor(x,y)
    #py.show()

    for i in range(0,len(l)):
        P = PSF_pix(x,y)[0]
        S = source_flux()[0][i]*t*P
        B = np.ones((len(xx),len(yy)))*BG()[3][i]*t
        Var = S + B + inst.dI[i]*t + inst.RON[i]**2
        #print np.sum(np.sqrt(Var))
        
        w = 1/Var
        norm = np.sum(w*P**2.)
        #w = (P**2./Var)/norm
        #print w[12,13]
        #print w[5,5]
        #print w[1,2]
        F = np.sum(w*S*P)/norm
        #print F
        #print source_flux(inst,color)[0][i]*t
        noise = np.sum(np.sqrt(Var/P**2))
        #print noise
        SR = F/noise
        SNR.append(SR)
    
    return SNR 
'''
