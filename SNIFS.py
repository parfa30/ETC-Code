#!/usr/bin/env python
from numpy import interp
from scipy import interpolate

from ETC_inputs import *
from Instrument import * 

class SNIFS(Instrument):
    def __init__(self):
        Instrument.__init__(self,'/Users/parkerf/Documents/Research/ETC/ETC/inst_parfiles/SNIFS.txt')
        
        self.option = 'dich_b' #This option refers which dichroic we think is used in SNIFS. There was some question as to whether it was dich_b or dich_a 
        self.b_channel = np.arange(self.par['low'],self.par['middle'],self.par['dlambda'])
        self.r_channel = np.arange(self.par['middle'],self.par['high'],self.par['dlambda'])

        self.res = np.append(self.par['b_res']*np.ones(len(self.b_channel)),self.par['r_res']*np.ones(len(self.r_channel)))
        self.dI = np.append(self.par['b_dI']*np.ones(len(self.b_channel)), self.par['r_dI']*np.ones(len(self.r_channel)))
        self.RON = np.append(self.par['b_RON']*np.ones(len(self.b_channel)), self.par['r_RON']*np.ones(len(self.r_channel)))
        
    def grid(self):
        #This helps define the MLA so that we can do PSF photometry
        x_num = self.par['xx']*self.par['spax_w'] 
        self.x_num = np.arange(0,x_num,self.par['spax_w'])
        self.x0 = x_num/2.
        y_num = self.par['yy']*self.par['spax_w']
        self.y_num = np.arange(0,y_num,self.par['spax_w'])
        self.y0 = y_num/2.
        return self.x_num,self.x0,self.y_num,self.y0
        
    def n(self,area):
        #this needs ot be thought through. THere is no I.area
        return (area/self.spax())/self.pix #number of pixels covered by PSF
    
    def Transmission(self,color):
        if color == 'b':
            self.path_ = '/Users/parkerf/Documents/Research/ETC/transmissions/dich_b/Blue_'
        elif color == 'r':
            self.path_ =  '/Users/parkerf/Documents/Research/ETC/transmissions/dich_b/Red_'
        
        Color = os.listdir(self.path_+self.option) 
        Ttotal = 1
        Tlist = []
        
        #This part used if grism model is judged inaccurate
        #if color == 'blue':
        #    Ttotal *= blue_grism(l)*0.8
        #    Tlist.append(blue_grism(l)*0.8)8
        #elif color == 'red':
        #    pass
        
        for element in Color: 
            if element == '.DS_Store':  
                pass
            elif element == 'blue_grating.txt':
                data = np.loadtxt(self.path_+self.option+'/'+element)
                ind1 = np.where( data[:,0] >= 4050 ) 
                ind2 = np.where( data[:,0] <4050 )
                data[:,1][ind1] *= 0.8 #Since the grating model isn't understood perfectly, This is used to adjust the tput manually until it looks like we think it should
                data[:,1][ind2] *= 0.7 #The data is split to deal with a particular spike in the grism 
                self.percent(data)
                output = interpolate.interp1d(data[:,0],data[:,1])
                Ttotal *= ((output(wave))) #np.ones(len(output(l))) if you want to take out the dependence on the grating completely
                Tlist.append((output(wave))) #np.zeros(len(output(l)))
               
            elif element == 'blue_QE.txt':
                data = np.loadtxt(self.path_+self.option+'/'+element) #for some reason this file seems to act up
                self.percent(data)
                output = interpolate.interp1d(data[:,0],data[:,1])
                Ttotal *= output(wave)
                Tlist.append(output(wave))
            
            else:
                data = np.loadtxt(self.path_+self.option+'/'+element)
                self.percent(data)
                output = interpolate.interp1d(data[:,0],data[:,1])
                Ttotal *= output(wave)
                Tlist.append(output(wave))
                  
        return Ttotal, Tlist 
    
    def Tput(self):
        Blue = self.Transmission('b')[0]
        Red = self.Transmission('r')[0]
        return Blue+Red
        
    def plot_Tput_components(self,color):
        element_list_b=self.Transmission('b')[1]
        element_list_r=self.Transmission('r')[1]
        
        from itertools import cycle
        lines = ["-","--","-.",":","-"]
        linecycler = cycle(lines)
        fig = py.figure(figsize=[7,5])
        ax  = fig.add_axes([0.15,0.1,0.8,0.8])
        py.title('Transmission (SNIFS+Tel) w/'+self.option)
        py.ylabel('Throughput %')
        py.xlabel('Wavelength ($\AA$)')
        py.ylim(0,1.0)
        #py.xlim(3000,10000)
        if color.lower() == 'b':
            for element in element_list_b:
                py.plot(wave,element,next(linecycler), color = 'b')
        elif color.lower() == 'r':
            for element in element_list_r:
                py.plot(wave,element,next(linecycler), color = 'r')
        else:
            print "Choose a color (b/r)"  
            
        py.show() 
    
    def plot_Tput(self):
        BTput = self.Transmission('b')[0]
        RTput = self.Transmission('r')[0]
        T = self.tel_tput()
        fig = py.figure(figsize=[7,5])
        ax  = fig.add_axes([0.15,0.1,0.8,0.8]) 
        l1,l2 = ax.plot(wave,BTput*T,'b',wave,RTput*T,'r') 
        py.title('Transmission of '+self.name+' + Tel')
        py.ylabel('Throughput %')
        py.xlabel('Wavelength ($\AA$)')
        py.ylim(0,0.5)
        py.grid()
        py.minorticks_on()
        
