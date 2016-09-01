# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 15:00:27 2016

@author: claire
"""

"""Main HVSR on seismes """



from Events_caracterisation import *
from Events_Selections import *
from Signals_Analysis import *
from Signals_Pre_Processing import *
from Plot_results import *
from Plot_Waveforms import *


jJultoplot = '110'
LStations = ['1','2','3','4','5','6','7']
#LStations = ['1']
MlMin= 3 
MlMax = 9
DistMin = 0
DistMax = 200
File = '/home/rault/PHD/Data/Event_stations_Caracteristics.txt'


Lenstation = np.zeros(8)
for Station in LStations : 
    #select ligne where stations 
    EventCaractSelec, shape = SelectEventMLDist(File, Station, MlMin, MlMax, DistMin, DistMax)
    LHVSR = []
    Lfr = []
    outfile = '/home/rault/PHD/Results/HVSR/HVSR_plots/HVSR_Seismes%s.eps'%Station
    outfileV = '/home/rault/PHD/HVSR/HVSR_plots/HVSR_Seismes%sv.eps'%Station
    outfiletxt = '/home/rault/PHD/Results/HVSR/HVSR_Arrays/HVSR_Seismes%s.npy'%Station
    outfiletxtfr = '/home/rault/PHD/Results/HVSR/HVSR_Arrays/HVSR_frh_Seismes%s.npy'%Station
    print len(EventCaractSelec), 'event caract'
    A = len(EventCaractSelec)
    for j in xrange(len(EventCaractSelec)):
        
        #1.0 take event characteristics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       Station, Year,jJul, Hour, Second = ConvertDatestr(EventCaractSelec[j])
      
       
       # 2.0 read file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       print jJul,  'read'      
       if not  os.path.exists('/home/burtin/DATA/LinTianShan/Seismic_Data/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHN.SAC'%(Year,jJul,Station,Year,jJul,str(Hour))) : 
           print 'file not existing' , jJul, Hour 
           A = A-1
           continue
       
       if jJul == jJultoplot: 
           plot =True
       else : 
           plot= False
       st=Read_event(Year,jJul,Hour,Second, Station, plot)
       
       #3.0 instrumental correction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       stcorrect = Stream_Correction(st, Station, False)
       del st
     
       
       #4.0 Trim signal between p theoirical arrival and end of signal defined as LTA/LTA ratio == 0.5~~~~~~~~~~
       st = Trim(stcorrect,Second,2, 0.5,10,40,False,plot)
       del stcorrect

       # 5.0 Vertical and Horizontal component~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       trH,trZabs = RealtrH(st,plot)
  
       
       #6.0 Smooth Spectrum: bandwidth = 110 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     
       frv, SV = SmoothSpectrum(trZabs, plot)
  
       frh, SH = SmoothSpectrum(trH, plot)
       del trZabs
       del trH
     
       
       #5.0 HSVR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
       Hvsr = HVSR(SH,SV)
       LHVSR.append(Hvsr)
       Lfr.append(np.array(frh))
                 
    AHVSR = np.asarray(LHVSR)
    Afr = np.asarray(Lfr)
    # 7.0 mean and standard variation of HVSR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print 'Mean'
    
    MeanSTDHVSR(AHVSR,Afr,outfile)

    np.save(outfiletxt,AHVSR)
    np.save(outfiletxtfr,frh)
    #record the number of event taken per stations
    Lenstation[int(Station)] = A
    
print 'number of event per station ', Lenstation
