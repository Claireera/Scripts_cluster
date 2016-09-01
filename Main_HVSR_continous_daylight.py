# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 15:15:27 2016

@author: claire
"""

"""Main HVSR continus 
~~~~~~~~~~~~~~~~~~~~~

HVSR is computed continouly during a given period 
HVSR is the mean of the hourly HVSR calculated along the day for given hours ex [0,5][22,23]
"""
from Events_caracterisation import *
from Events_Selections import *
from Signals_Analysis import *
from Signals_Pre_Processing import *
from Plot_results import *
from Plot_Waveforms import *
Year ='15'

#Lmonth = [range(90,119),range(120,150),range(151,180),range(181,211),range(212,242),range(243,272),range(273,303),range(304,333),range(305,364)]
Lmonth = range(59,364)
LHour = range(0,5)+range(14,24) #day light
plageH = '5_20'
#LHour = range(14,22) #night
#plageH = '22_5'

jJultoplot=False
Hourtoplot = False
#LStations = ['1','2','3','4','5','6','7']
LStations = ['7']

for Station in LStations:
    m=0

    AllHVSR = []
    Allfrh = []
    AllJJul = []
    m+=1
    for jJul in Lmonth:
	    	LHVSR = []
	    	Lfr = []
	    	outfile = '/home/rault/PHD/Results/HVSR/HVSR_plots/HVSR_Continuous_%s_H_%s.eps'%(Station,plageH )
		outfiletxt = '/home/rault/PHD/Results/HVSR/HVSR_Arrays/HVSR_Continuous_HVSR_%s_H_%s.npy'%(Station,plageH)
		outfilefrtxt = '/home/rault/PHD/Results/HVSR/HVSR_Arrays//HVSR_Continuous_fr_%s_H_%s.npy'%(Station,plageH)
		outfilejUltxt = '/home/rault/PHD/Results/HVSR/HVSR_Arrays//HVSR_Continuous_jJul_%s_H_%s.npy'%(Station,plageH)
		

		for Hour in LHour:
			if Hour<10: 
				Hour='0'+str(Hour)
		    	else : 
				Hour=str(Hour)
			if jJul<100:
				jJul = '0'+str(jJul)
		    	# 1.0 read file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		    	if not os.path.exists('/home/burtin/DATA/LinTianShan/Seismic_Data/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHN.SAC'%(Year,str(jJul),Station,Year,str(jJul),Hour)) : 
			       	print 'file not existing, Jjul :  ' , jJul, 'Hour', Hour,'Station', Station
			       	Hvsr = np.empty((180299,0,))
			       	Hvsr[:]=np.NAN
				jJul = int(jJul)
			       	#frh = Hvsr
		       		continue
			if int(jJul)==101 and int(Station)==7:
				print 'jjul 101'
				continue
			if int(jJul)>244 and int(jJul)<376 and int(Station)==7:
				print 'jul 101'
				continue
			print 'Julian day ', jJul, Hour
		   	if jJul == jJultoplot and Hour==Hourtoplot: 
				plot =True
		    	else : 
				plot= False
		    
		    	st=Read_event(Year,str(jJul),Hour,3599, Station, plot)
		   	jJul = int(jJul)
		   	# 2.0 instrumental correction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		    	stcorrect = Stream_Correction(st, Station, plot)
		    	del st
		    	st = stcorrect.copy()
		    	del stcorrect
		    	# 3.0 Vertical and Horizontal component~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		    	
		    	trH,trZabs = RealtrH(st,plot)
		   
		    	# 4.0 Smooth Spectrum: flat :300 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		        
		    	frv, SV = SmoothSpectrum(trZabs, plot)
		    	frh, SH = SmoothSpectrum(trH, plot)
		 
		    	del trZabs
		    	del trH
		       
		    	# 5.0 HSVR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	          
		    	Hvsr = HVSR(SH,SV)
		    	LHVSR.append(Hvsr)
		    	Lfr.append(frh)
   		AHVSR= np.array(LHVSR)
    		Afr = np.array(Lfr)
   		Mean = np.nanmean(AHVSR,axis=0).tolist()
   		fr = np.nanmean(Afr, axis=0).tolist()
   		if isinstance(Mean, float): 
        		Mean = np.zeros(50299).tolist()
    		if isinstance(fr, float): 
        		print 'nodat'
        		print fr, "ll"
    		else:
        		Allfrh.append(fr)
        		Fr=fr
        	AllHVSR.append(Mean)
        	AllJJul.append(jJul)
    # 6.0 Save HVSR means fr and Julain day arrays for a given station ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AllHVSR = np.asarray(AllHVSR).T   
    np.save(outfiletxt,AllHVSR)
   
    PlotSpectrogram(Fr, AllHVSR, AllJJul, outfile,colormap="jet")


