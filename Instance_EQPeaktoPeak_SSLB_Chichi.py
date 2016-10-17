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

class PtPearthquake : 
   
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
        self.PeaktoPeak = 0
        self.PGA = 0
        self.SNR = 0
        self.valid = 0
        self.PGV = 0
        self.Arias =0 
        self.second_P = 0
        self.second_S = 0
        self.Component = component
        self.trFilt =[]
        self.trTrim = []
        self.Enveloppe =0
      
        return
    
    def __str__(self):
        return "St : {} - ML : {} - Rdist:  {} - Time :{}" .format(self.Station,self.magnitude,self.depth, self.Rdistance,self.time)
        
    def Save(self):
        """save class as self.name.txt"""
        file = open('/home/rault/PHD/Results/Peak_to_Peak/PeaktoPeak_Data/'+str(self.time[0])+'/'+'St'+str(self.Station)+'/St'+str(self.Station)+str(self.Component)+'/'+self.nameEQ+'_'+self.Station+'_'+self.Component+str(self.frband[1])+'.txt','w')
        file.write(cPickle.dumps(self.__dict__))
        file.close()
        print 'SAVED',str(self.time[0])+'/'+'St'+str(self.Station)+'/St'+str(self.Station)+str(self.Component)+'/'+self.nameEQ+'_'+self.Station+'_'+self.Component+str(self.frband[1])
        return 
        
                
    def Dict(self, dataPickle):
        self.__dict__ = cPickle.loads(dataPickle)
        return 
        
            
        
    def FStaturation(self):
        
        """define if the signal of a given component of an EQ is saturated or not"""
        self.Satured= Saturation(self.tr)
        return

    def FPGA(self):
        """calculate the maximum of acceleration of each component of a given EQ"""
   
        self.PGA= np.max(np.abs(np.gradient(self.tr.data)))
        return
        
    def FPGV(self):
        """calculate the maximum of velocity of each component of a given EQ"""
        self.PGV= np.max(np.abs(self.tr.data))
        return
    def FArias(self):
        """calculate the Arias of each component of a given EQ"""
       
        self.Arias= self.tr.stats.sampling_rate*simps(np.gradient(self.tr.data)**2)*pi/2*9.81
        return
    def FEnveloppe(self):
        """calculate the maximum of the eneveloppe of the signal calculated nvelope is determined by adding the squared amplitudes of the function and it’s Hilbert-Transform and then taking the square-root. (See [Kanasewich1981]))"""
        
        self.Enveloppe = np.max(obspy.signal.filter.envelope(self.tr.data))
        return
        
    def Filter(self):
        """filter the signal with a butterworth between for a given frband for a frequency band of 2Hz with moving of 1Hz """
        trfilt = self.tr.copy()   
        Fmin =  self.frband[0]
        Fmax =  self.frband[1]
        print 'filt'
        self.trFilt = trfilt.filter('bandpass',freqmin=Fmin,freqmax=Fmax,corners=4,zerophase=True) 
        print "filt the "
        return      
        
    def Filter_gaussian(self,std):
        """filter the signal with a gaussian between for a given central frequency (mu) and a std """
        trfilt = self.tr.copy()   
        std = std
        print 'filtgaussian'
        self.trFilt = tr_Gaussian_filter2(self.tr, std, self.mu)
        print "filtg aussain "
        return      
        
        
    def SNR(self,LimSNR):
        """return SRN calculated as baillard and if the sigal is valid ie SRN above minSNR
        input
        limSRN: SNR limite above the signal is valid
        output
        self.SNR 
        self.valid
        """
        self.valid, self.SNR= SNR(self.tr, self.trTrim, LimSnR, self.second_P) 
        
        return 
        
        
    def PtP_envelope(self):
        
         """return a list of array type ([Maxpeak],[PeaktoPeak],[Sumenvelope]) at each turn i function made to share work between the threads
         input : 
         output :instance attributes are defined : Envelope_Itegral, and PtP for each component of the stream 
         exemple : EQ.PtP_envelope"""
    
         #2.3. Peak to Peak ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         self.PeaktoPeak, ptpFirstR = PeaktoPeak(self.tr,1,False)
  
         return