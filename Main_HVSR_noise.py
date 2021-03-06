# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 19:21:16 2016

@author: claire
"""

"""Main HVSR noise"""
from Events_caracterisation import *
from Events_Selections import *
from Signals_Analysis import *
from Signals_Pre_Processing import *
from Plot_results import *
from Plot_Waveforms import *

LjJday =[215,240]

 
LStations = ['1','2','3','4','5','6','7']
Year = '15'
jJultoplot =258
Hourtoplot = '15'

for Station in LStations : 
       
    for jJul in LjJday:
        LHVSR = []
        Lfr = []
        outfile = '/home/rault/PHD/Results/HVSR/HVSR_plots/HVSR_Noise_jJul_%s_%s.eps'%(str(jJul),Station)
        outfiletxt = '/home/rault/PHD/Results/HVSR/HVSR_Arrays/HVSR_Noise_jJul_%s_%s.npy'%(str(jJul),Station)
        outfilefrtxt = '/home/rault/PHD/Results/HVSR/HVSR_Arrays/HVSR_Noise_jJul_fr%s_%s.npy'%(str(jJul),Station)
   	if jJul==215 :
	    Lhour = np.concatenate((range(0,9),range(10,13),range(14,24)))
	elif jJul==240:
	    Lhour =np.concatenate((range(0,5),range(6,7),range(8,24)))
	for Hour in xrange(0,24,):
            if Hour<10: 
                Hour='0'+str(Hour)
            else : 
                Hour=str(Hour)
            # 1.0 read file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if not  os.path.exists('/home/burtin/DATA/LinTianShan/Seismic_Data/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHN.SAC'%(Year,str(jJul),Station,Year,str(jJul),Hour)) : 
               print 'file not existing' , jJul, Hour 
               continue
	    if int(jJul)==101 and int(Station) ==7:
	       continue
            if jJul == jJultoplot and Hour==Hourtoplot: 
                plot =True
            else : 
                plot= False
            
            st=Read_event(Year,str(jJul),Hour,3599, Station, plot)
           
            # 2.0 instrumental correction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            stcorrect = Stream_Correction(st, Station, plot)
            del st
            st = stcorrect.copy()
            del stcorrect
            # 3.0 Vertical and Horizontal component~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            trH,trZabs = RealtrH(st,plot)
           
            # 4.0 Smooth Spectrum: bandwidth = 110 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            print 'SMOOTH', jJul, '_',  Hour       
            #KonnoSmoothSpec(trZ,110,plot)
            frh, SH =SmoothSpectrum(trH, plot)
            #KonnoSmoothSpec(trH,110,plot)
	    frZ, SV =SmoothSpectrum(trZabs, plot)

            del trZabs
            del trH
               
            # 5.0 HSVR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            print 'HVSR'            
            Hvsr = HVSR(SH, SV)
            LHVSR.append(Hvsr)
            Lfr.append(frh)
            
        ALHVSR = np.asarray(LHVSR)   
        Afrh=np.asarray(Lfr)
        # 6.0 mean and standard variation of HVSR  plot and save plot
        MeanSTDHVSR(ALHVSR,Afrh, outfile)
        np.save(outfiletxt,ALHVSR)
  
         
