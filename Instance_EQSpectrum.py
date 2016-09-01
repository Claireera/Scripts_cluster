# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 17:37:35 2016

@author: claire
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 16:00:25 2016

@author: claire
"""

"""Build Earthquake Peak to Peak instance """
import cPickle
from Events_Selections import *
from Signals_Analysis import *
from Signals_Pre_Processing import *
from Plot_results import *
from Plot_Waveforms import *
import numpy as np
from scipy.integrate import simps
from Signal_Wave_Picking import *


class Spectrumearthquake : 
   
    """Ptpearthquake is an instance that give 
    EQ1= PtPearthquake(name,magnitude,depth,Rdistance, Hdistance,Lat,Long,Az,time)
    PtP1 = EQ1.N"""
    
    def __init__(self,name,magnitude,depth,Rdistance, Lat,Long,Az,time,Station,tr, component):
        self.nameEQ = name
        self.depth = depth
        self.magnitude =magnitude
        self.Rdistance = Rdistance
        self.Coord = (Lat,Long)
        self.Az = Az
        self.time = time
        self.Station = Station
        self.tr=tr
        self.Satured=0
        self.spectrum = []
        self.trFilt =[]
        self.trTrim = []
        self.frequencies = []
        self.component = component
        self.second_P = 0
        self.second_S = 0
    
    def __str__(self):
        return "St : {} - ML : {} - Rdist:  {} - Time :{}" .format(self.Station,self.magnitude,self.depth, self.Rdistance,self.time)
        
    def Save(self):
        """save class as self.name.txt"""
        file = open('/home/rault/PHD/Results/Spectrum/'+self.nameEQ+'_'+self.Station+'_'+self.component+'.txt','w')
        file.write(cPickle.dumps(self.__dict__))
        file.close()
        return 
        
                
    def Dict(self, dataPickle):
        self.__dict__ = cPickle.loads(dataPickle)
        return 
        
        
    def FStaturation(self):
        
        """define if the signal of a given component of an EQ is saturated or not"""
        self.Satured= Saturation(self.tr)
                
    def Filter(self):
        """filter the signal with a butterworth between for a given frband for a frequency band of 2Hz with moving of 1Hz """
        trfilt = self.tr.copy()   
        Fmin = 0.01 + self.frband[0]
        Fmax = 2+ self.frband[1]
        print 'filt'
        self.trFilt = trfilt.filter('bandpass',freqmin=Fmin,freqmax=Fmax,corners=4,zerophase=True) 
        print "filt the "
        return      
        
    def SNR(self,LimSNR):
        """return SRN calculated as baillard and if the sigal is valid ie SRN above minSNR
        input
            - limSRN: SNR limite above the signal is valid
        output
            - self.SNR 
            - self.valid
        """
        self.valid, self.SNR= SNR(self.tr, self.trTrim, LimSnR, self.second_P) 
        
        return 
        
    def Spectrum(self):
        
         """return smooth spectrum of the EQ trace for the component considered"""
              
         self.frequencies, self.Spectrum = SmoothSpectrum(self.trTrim, False)

         return


      
