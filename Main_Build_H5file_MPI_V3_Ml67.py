# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 14:15:33 2016

@author: claire
"""
import h5py

# Initialize the communicator
import mpi4py
from mpi4py import MPI

Comm = MPI.COMM_WORLD

me = Comm.Get_rank()  # The process ID (integer 0-3 for 4-process run)
size = Comm.Get_size()  # The number of processes

timestart=MPI.Wtime()

#t2 = time.clock() 

#from Events_caracterisation import *
from Events_Selections import *
from Signals_Pre_Processing import *
from Plot_results import *
from Plot_Waveforms import *
from Instance_EQPeaktoPeak import *
from Signal_Wave_Picking import *
import cPickle
import json
import calendar
"""0. Parameters Initialisation"""

# limit SNR above signal is valid

#Magnitude,min and max, and Distance min and max
MlMin= 6
MlMax = 7
DistMin = 0
DistMax = 200
jJultoplot = 206
refStation = '5'
LStations = ['1','2','3','4','5','6','7']
#component 
LComponents = ['N','R','E','T','Z','H']
# 0. limit SNR above signal is valid
LimSNR = 2
"""1. Open Earthquakes file comtaining eqrthqukes fetures station, component and  Selection (in a separate script) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""

with open('/home/rault/PHD/Data/List_EQ_St_Comp_Freq_Ml_%s_%s_distmax_%s_05.txt'%(MlMin,MlMax,DistMax)) as f:
    LEqStCompFreq =json.load(f)
f = h5py.File('/home/rault/PHD/Results/Parallel_PeaktoPeak_test2_V3_Ml_%s_%s_distmax_%s_05.hdf5'%(MlMin,MlMax,DistMax), 'w',driver= 'mpio', comm=Comm,libver='latest')
"""5. Build database : h5 file structure~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""
for EqSta in LEqStCompFreq:
    # Get event characteristics
    EQname, Station, Year,jJul, Hour, Secondp,Seconds, Ml,Depth, Rdistance, Lat,Long, Az,BAz,component,number, freqband = EqSta
  
    if EQname not in f : 
        grpEQ = f.create_group(EQname)
        #print 'create group ', EQname
        # Save the metadata (lon, lat, etc) of the earthquake ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #longitude
        grpEQ.attrs['Long'] = Long
        #latitude
        grpEQ.attrs['Lat'] = Lat
        #Magnitude
        grpEQ.attrs['Ml'] = Ml
        #depth
        grpEQ.attrs['Depth'] = Depth
        #Time
        grpEQ.attrs['Year'] = Year 
        grpEQ.attrs['JJul'] = jJul 
        grpEQ.attrs['Hour'] = Hour 
        grpEQ.attrs['Secondp'] = Secondp  
        grpEQ.attrs['Seconds'] = Seconds  
        #create a groupe for station (each line correspond to a EQ station couple)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        grpSt = f[EQname].create_group('St{}'.format(Station))
        #save the metadata of the Station 
        grpSt.attrs['RDist'] = Rdistance
        grpSt.attrs['Az'] = Az
        grpSt.attrs['BAz'] = BAz
        for Component in LComponents:
            if 'PtP_{}'.format(Component) not in f[EQname]['St{}'.format(Station)]:
                dStPtP = f[EQname]['St{}'.format(Station)].create_dataset('PtP_{}'.format(Component),(12,), dtype='f')
                dSecondS = f[EQname]['St{}'.format(Station)].create_dataset('SecondS_{}'.format(Component),(12,), dtype='f')
                dSecondP = f[EQname]['St{}'.format(Station)].create_dataset('SecondP_{}'.format(Component),(12,), dtype='f')
                dfreq = f[EQname]['St{}'.format(Station)].create_dataset('Freq_{}'.format(Component),(12,), dtype='f')
                dStPGA = f[EQname]['St{}'.format(Station)].create_dataset('PGA_{}'.format(Component),shape=(1,),dtype='f')
                dStPGV = f[EQname]['St{}'.format(Station)].create_dataset('PGV_{}'.format(Component),shape=(1,),dtype='f')
                dStArias = f[EQname]['St{}'.format(Station)].create_dataset('Arias_{}'.format(Component),shape=(1,),dtype='f')
                dStSatured = f[EQname]['St{}'.format(Station)].create_dataset('Satured_{}'.format(Component),shape=(1,),data=True)
                dStvalid = f[EQname]['St{}'.format(Station)].create_dataset('valid_{}'.format(Component),shape=(1,),data=False)
                dStSNR = f[EQname]['St{}'.format(Station)].create_dataset('SNR_{}'.format(Component),shape=(1,),dtype='f')
                dStart= f[EQname]['St{}'.format(Station)].create_dataset('Start_time_{}'.format(Component),shape=(1,),dtype='f') 
                dSend= f[EQname]['St{}'.format(Station)].create_dataset('End_time_{}'.format(Component),shape=(1,),dtype='f') 
                dStexist = f[EQname]['St{}'.format(Station)].create_dataset('Exist_{}'.format(Component),shape=(1,),data=True)
                dStProcessor = f[EQname]['St{}'.format(Station)].create_dataset('Processor_{}'.format(Component),(12,), dtype='f')
                dStrank = f[EQname]['St{}'.format(Station)].create_dataset('rank_{}'.format(Component),shape=(12,),dtype='f')
                
    if 'St{}'.format(Station) not in f[EQname]: 
            #print 'creation group Station ' , Station , 'EQ ', EQname
            grpSt = f[EQname].create_group('St{}'.format(Station))
            #save the metadata of the Station 
            grpSt.attrs['RDist'] = Rdistance
            grpSt.attrs['Az'] = Az
            grpSt.attrs['BAz'] = BAz
            
            for Component in LComponents:
                if 'PtP_{}'.format(Component) not in f[EQname]['St{}'.format(Station)]:
                    dStPtP = f[EQname]['St{}'.format(Station)].create_dataset('PtP_{}'.format(Component),(12,), dtype='f')
                    dSecondS = f[EQname]['St{}'.format(Station)].create_dataset('SecondS_{}'.format(Component),(12,), dtype='f')
                    dSecondP = f[EQname]['St{}'.format(Station)].create_dataset('SecondP_{}'.format(Component),(12,), dtype='f')
                    dfreq = f[EQname]['St{}'.format(Station)].create_dataset('Freq_{}'.format(Component),(12,), dtype='f')
                    dStPGA = f[EQname]['St{}'.format(Station)].create_dataset('PGA_{}'.format(Component),shape=(1,),dtype='f')
                    dStPGV = f[EQname]['St{}'.format(Station)].create_dataset('PGV_{}'.format(Component),shape=(1,),dtype='f')
                    dStArias = f[EQname]['St{}'.format(Station)].create_dataset('Arias_{}'.format(Component),shape=(1,),dtype='f')
                    dStSatured = f[EQname]['St{}'.format(Station)].create_dataset('Satured_{}'.format(Component),shape=(1,),data=True)
                    dStvalid = f[EQname]['St{}'.format(Station)].create_dataset('valid_{}'.format(Component),shape=(1,),data=False)
                    dStSNR = f[EQname]['St{}'.format(Station)].create_dataset('SNR_{}'.format(Component),shape=(1,),dtype='f')
                    dStart= f[EQname]['St{}'.format(Station)].create_dataset('Start_time_{}'.format(Component),shape=(1,),dtype='f') 
                    dSend= f[EQname]['St{}'.format(Station)].create_dataset('End_time_{}'.format(Component),shape=(1,),dtype='f') 
                    dStexist = f[EQname]['St{}'.format(Station)].create_dataset('Exist_{}'.format(Component),shape=(1,),data=True)
                    dStProcessor = f[EQname]['St{}'.format(Station)].create_dataset('Processor_{}'.format(Component),(12,), dtype='f')
                    dStrank = f[EQname]['St{}'.format(Station)].create_dataset('rank_{}'.format(Component),shape=(12,),dtype='f')
print 'database built'
print 'fclosing'          
# Close h5file
f.close()
# Finalize
MPI.Finalize()
# End time
timeend = MPI.Wtime()
print ('---------------------------------------------/n','time start : ',timestart,'/n','time end : ',timeend,'\n---------------------------------------------/n')
