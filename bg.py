#!/usr/bin/env python
import os
import glob
import sys
from numpy import interp
import scipy.io as sio

from ETC_inputs import *
from Instrument import *



# All atm lines are for airmass = 1. Will change later
# Most of these files are from Gemini - except for the Atm_ext in the optical (Buton paper) - In IR, water vapor column depth = 1.6mm as a standard. Would like to be able to change this later on

wave = np.arange(I.lower,I.upper,I.step)
l = wave

def Ext_Atm():

    #Atm. extinction from Buton paper. Covers 3000A - 10300A
    optical = np.loadtxt('/Users/parkerf/Documents/Research/ETC/transmissions/med_atm_ext_mk.txt')
    Oext = interp(wave,optical[:,0],optical[:,1])
    
    #The following are mean extinction values for mosfire bands up to ~24000A. Taken from Gemini website. They do not use these numbers in their ITCs.
    J = np.arange(10301,13519,1)
    Jext = interp(wave,J,0.015*np.ones(len(J)),left=0,right=0)
    H = np.arange(13520,19539,1)
    Hext = interp(wave,H,0.015*np.ones(len(H)),left=0,right=0)
    K = np.arange(19540,23970,1)
    Kext = interp(wave,K,0.015*np.ones(len(K)),left=0,right=0)

    bands = (Oext,Jext,Hext,Kext)
    Text = []
    
    for b in bands:
        Text.append(10**(-.4*I.amass*b))

    Ext = np.prod(Text,axis=0)
        
    return Ext  

def Atm_trans_lines():
    #These are the IR transmission lines. 0.9-5.6 Microns. From Gemini Website
    wvap = str(int(float(I.wvap)*10))
    amass = str(int(float(I.amass)*10)) 
    data = np.genfromtxt('http://www.gemini.edu/sciops/ObsProcess/obsConstraints/atm-models/'+I.loc+'trans_zm_'+wvap+'_'+amass+'.dat')
    Trans = interp(l,data[:,0]*10000,data[:,1],left=1,right=1)
    return Trans

def Atm_emiss_lines():

    #IR lines from Gemini webiste. 0.9 - 5.6um. photons/s/arcsec^2/A/m^2
    wvap = str(int(float(I.wvap)*10))
    amass = str(int(float(I.amass)*10))
    ir_data = np.genfromtxt('http://www.gemini.edu/sciops/ObsProcess/obsConstraints/atm-models/'+I.loc+'_skybg_zm_'+wvap+'_'+amass+'_ph.dat')
    IR_data = interp(l,ir_data[:,0]*10,ir_data[:,1],left=0,right=0)

    IR = IR_data*inst.apt()*(0.0001)*inst.spax() #gives photons/sec/arcsec^2/A
    
    #Optical lines. Need to be very wary of the area 3740-3760. There could be double counting in this area. Seems low enough not to worry about though. 10^-16^erg/s/cm2/A
    o1='/Users/parkerf/Documents/Research/Hubble/J_A+A_407_1157/table4.dat' #Sky emission lines from 3140 to  3760 {AA}
    o2='/Users/parkerf/Documents/Research/Hubble/J_A+A_407_1157/table5.dat' #Sky emission lines from 3740 to  4850 {AA}
    o3='/Users/parkerf/Documents/Research/Hubble/J_A+A_407_1157/table6.dat' #Sky emission lines from 4850 to  5770 {AA}
    o4='/Users/parkerf/Documents/Research/Hubble/J_A+A_407_1157/table7.dat' #Sky emission lines from 5830 to  6700 {AA}
    o5='/Users/parkerf/Documents/Research/Hubble/J_A+A_407_1157/table8.dat' #Sky emission lines from 6700 to  8540 {AA}
    o6='/Users/parkerf/Documents/Research/Hubble/J_A+A_407_1157/table9.dat' #Sky emission lines from 8600 to 10430 {AA}

    files = (o1,o2,o3,o4,o5,o6)
    filenames = ['o1','o2','o3','o4','o5','o6']
    
    #opt = []
    #for f in files:
        #data = np.genfromtxt(f)
        #opt.append(interp(l,data[:,1],data[:,3],left=0,right=0))

    #Opt = np.sum(opt,axis=0)
    #OPT = (Opt*inst.apt()*(l/(h*c))*10**-16*inst.spax()) #gives photons/sec/arcsec^2/A
    
    opt = np.loadtxt('/Users/parkerf/Documents/Research/ETC/transmissions/skybg_50_10.txt')
    Opt = interp(l,opt[:,0]*10,opt[:,1],left=0,right=0)
    OPT = Opt*inst.apt()*(.0001)*inst.spax() #photons/s/arcsec^2/A
    #ind = np.where(l>9000.)
    #OPT[ind] = 0
    
    return OPT + IR

def optical_sky():
    #Not using this currently
    U = [3650., 660.]
    B = [4450., 940.]
    V = [5510., 880.]
    R = [6580., 1380.]
    bands = (U,B,V,R)
    
    bandwidth = []  
    for b in bands:
        bandwidth.append(np.arange(b[0]-(b[1]/2.), b[0]+(b[1]/2.), inst.res))
    
    if I.sky == 'darkest':
        mag = [21.3, 22.1, 21.3, 20.4]
    elif I.sky == 'dark':
        mag =  [20.7, 19.2, 20.9, 19.9]
    elif I.sky == 'grey':
        mag = [19.5, 17.3, 19.5, 19.1]
    elif I.sky == 'bright':
        mag = [18, 15,17.5, 17.9]

    Flux = []
    for i,bn in enumerate(bandwidth):
        flux = np.ones(len(bn))*10**(-0.4*mag[i]+48.6) #(c/l**2)*
        Flux.append(interp(l,bn,flux, left = 0, right = 0))
    
    return np.sum(Flux,axis=0)
    
       
                                                                                                            
def BG(): 
    #Output is photons/s/A/pixel. Galactic extinction?
    
    Atm = Atm_trans_lines()#Ext_Atm()*
    Tput = inst.Tput()#inst.tel_tput()*
    
    Therm = 10**(-17.)*((l/10000.)**(-4.))*((np.exp((h*c)/(k*inst.Ttel*l))-1.)**(-1.))*inst.therm_emiss()*(l/(h*c))*inst.spax()*Tput*inst.disp

    Zodi = 10**(-17.)*((l/10000.)**(-4.64))*((np.exp((h*c)/(k*5800*l))-1)**(-1))*inst.apt()*(l/(h*c))*inst.spax()*Tput*inst.disp#*Atm
    
    if inst.loc == 'g':
        Sky = Atm_emiss_lines()*Tput*inst.disp#Atm*
    elif inst.loc == 's':
        Sky = np.zeros(l)
    else:
        print "Where is the instrument located?"
    BGcomps = [Therm,Zodi,Sky]
    
    return np.sum(BGcomps,axis=0), BGcomps
    

def plot_BG():

    fig,ax1 = py.subplots(1,figsize=[7,5])
    l1,l2,l3,l4 = ax1.plot(l,BG()[0], l,BG()[1][0],'--',l,BG()[1][1],':',l,BG()[1][2])
    ax1.legend( (l1,l2,l3,l4),('Total','Thermal','Zodiacal','Sky'),loc=1)
    
    ax1.set_title('Sky Background (e/s/pixel): '+inst.name) 
    ax1.set_ylabel('Total')
    #ax3.set_ylabel('Red')
    py.minorticks_on()


