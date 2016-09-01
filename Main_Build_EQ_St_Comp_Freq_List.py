# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 11:36:56 2016

@author: claire

Build an array of earthquakes features, station, component, and frequencies calculated from 0.2 to 10."""  
from Events_Selections import *
from Signals_Analysis import *
from Signals_Pre_Processing import *
import json 
import numpy as np
import copy
#Magnitude,min and max, and Distance min and max
MlMin= 3
MlMax = 4
DistMin = 0
DistMax = 100
jJultoplot = 206
refStation = '5'
LStations = ['1','2','3','4','5','6','7']


# 1. Selection (in a separate script) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
File = '/home/rault/PHD/Data/Event_stations_Caracteristics_2015_2016Full2.txt' 
Outputfile = open('/home/rault/PHD/Data/List_EQ_St_Comp_Freq_Ml_%s_%s_distmax_%s.txt'%(MlMin,MlMax,DistMax), mode = 'w+')

# This selection returns a list of tuples (station_n, earthquake_m) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A_St_EQ = Array_St_EQ_MlDist(File, MlMin, MlMax, DistMin, DistMax)

Ldata = []
for line in A_St_EQ:
    #Earthquakes metadata (one line correspond to one EQ and one station)
    station, year,jJul, hour, Secondp,Seconds,ml,depth, Rdistance, Lat,Long, Az = ConvertDatestr(line)
    BAz = line[9]
    EQname = str(year)+'_'+str(jJul)+ '_'+str(hour)+'_'+str(int(Secondp/60))   
#    LEQ.append(EQname)
#    LSt.append(station)
    Ldata.append([EQname, station, year,jJul, hour, Secondp,Seconds,ml,depth, Rdistance, Lat,Long, Az,BAz])

print 'number of  earthquakes selected with Ml %s max distence is %s'%(MlMax,MlMax), A_St_EQ.shape
#frequencies 
frqmin = [k/10. for k in range(2,110,10)]
frqmax = range(2,13,1)
Lfrqtuple = np.column_stack([frqmin,frqmax]).tolist()

#component 
LComponents = ['N','R','E','T','Z','H']
#list of couple frq component
LfreqTuplescompo = []
for j in xrange(len(LComponents)):
    for i in xrange(len(Lfrqtuple)):
        LfreqTuplescompo.append([LComponents[j],i,Lfrqtuple[i]])
#build a list of earthquakes and its features, plus station and frequencies : thus one worker will work for a given EQ, station, component and frequency   
LEqStCompFreq = Ldata
LEqStCompFreq2= []
for EQSt in LEqStCompFreq :
    
    for Compofrq in LfreqTuplescompo:
        EQst_copy = copy.copy(EQSt)
        print EQSt
        for j in xrange(len(Compofrq)) :
            EQst_copy.append(Compofrq[j])  
            print EQst_copy,j
        LEqStCompFreq2.append(EQst_copy)
# This selection returns a list of tuples (staftion_n, earthquake_m)
# Write sequence of lines at the end of the file.ll
json.dump(LEqStCompFreq2,Outputfile)
Outputfile.close()
