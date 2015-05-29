#!/usr/bin/env python
import numpy as np
from numpy import interp
from scipy import interpolate


#constants
h=6.62607*10**(-27.)                 #erg*s
c=2.99*10**(18.)                     #A/s
k=1.3806*10**(-16.)                  #erg/K
erg = 6.2415*10**(-11.)              #eV


from ETC_inputs import *
#The class Instrument contains mostly the general telescope definitions. It calls on instrument specific information which are held in different programs
class Instrument(object):
    def __init__(self,paramfile):

        #Parameters are saved in a file by the instrument name
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
        
        self.name = self.par['name'] #LRIS, MOSFIRE, SNIFS, WFIRST, GMOS_S
        self.type = self.par['type'] #optica/ir
        self.loc = self.par['location'] #ground(g), space(s)
        self.Ttel = self.par['Ttel'] #temperature of telescope
        self.mdeg = self.par['mirror_degrade'] #mirror degradation
        self.nmirr = self.par['mirrors'] #number of mirrors
        

    def percent(self,array):
        array[:,1] /= 100.

    def apt(self):
        #obscuration either in terms of a diameter or percentage. Max diameter allowed is 5m. Answer given in cm^2
        if self.par['obsc'] > 1.:
            return np.pi*(0.5*(self.par['dia']-self.par['obsc']))**2.
        elif self.par['obsc'] < 1.:
            return np.pi*((0.5*self.par['obsc']**2.)*(1.-self.par['obsc']))
        else: 
            return "Obscuration in terms of diameter of obstruction (max 500cm) or percentage of primary."
    
    def therm_emiss(self):
        #effective emissivity * emission area to be multiplied by thermal blackbody. In e/spax/s
        emiss_prim = 1. - (1.-self.par['emissivity'])**(self.par['obsc'])
        emiss_black = 1.
        emiss_other = 1.- (1.-self.par['emissivity'])**(self.par['other_mirrors'])
        if self.par['obsc'] > 1:
            primary = (np.pi*((self.par['obsc']-self.par['obsc'])/2.)**2.) * emiss_prim/((1.-self.par['emissivity'])**(self.par['obsc']))
            obsc = np.pi*(self.par['obsc']/2.)**2 * emiss_black
            other = (np.pi*((self.par['obsc'])/2.)**2.) * emiss_other/(1-emiss_other) #assume fills full solid angle
            return (primary + obsc + other)
        if self.par['obsc'] < 1:
            primary = np.pi*(self.par['obsc']/2.)**2*(1-self.par['obsc']) * emiss_prim/((1.-self.par['emissivity'])**(self.par['obsc']))
            obsc = self.par['obsc'] * emiss_black
            other = (np.pi*((self.par['obsc'])/2.)**2.) * emiss_other/(1-emiss_other)
            return (primary + obsc + other)
            

'''
This is for an instrument that has a flat efficiency. Don't need it right now.
    else:
        if color == 'blue':
            b = np.ones(len(inst.bL))*inst.eff
            bl = np.append(b,np.zeros(len(inst.rL)))
            B=np.array(zip(l,bl))
            b_inst = interpolate.interp1d(l,bl)
            return b_inst(l),0
        elif color == 'red':
            r = np.zeros(len(inst.bL))
            rl = np.append(r,np.ones(len(inst.rL))*inst.eff)
            R=np.array(zip(l,rl))
            r_inst = interpolate.interp1d(R[:,0],R[:,1])
            return r_inst(l), 0
        else:
            print 'Which color are you looking for (red/blue)?'

'''
