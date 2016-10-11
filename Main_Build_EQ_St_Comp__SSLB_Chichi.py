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
MlMax =4
DistMin = 0
DistMax = 100
jJultoplot = 206
refStation = '5'
LStations = ['1','2','3','4','5','6','7']

 
# 1. Selection (in a separate script) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#File = '/home/rault/PHD/Data/Event_stations_Caracteristics_2015_2016Full2.txt' 
File = '/home/rault/PHD/Data/Event_stations_Caracteristics_SSLB_Chichi_afterShocks' #BATS
Outputfile = open('/home/rault/PHD/Data/List_EQ_St_Comp_Chichi_afterShocks.txt', mode = 'w+')

# This selection returns a list of tuples (station_n, earthquake_m) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#A_St_EQ = Array_St_EQ_MlDist(File, MlMin, MlMax, DistMin, DistMax)
A_St_EQ = Array_St_EQ_MlDist_BATS(File, MlMin, MlMax, DistMin, DistMax) # Bats
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

#component 
LComponents = ['N','R','E','T','Z','H']

#build a list of earthquakes and its features, plus station and frequencies : thus one worker will work for a given EQ, station, component and frequency   
LEqStComp = Ldata
LEqStComp2= []
for EQSt in LEqStComp :
    
    for Compo in LComponents :
        EQst_copy = copy.copy(EQSt)
        print EQSt
        
        EQst_copy.append(Compo)  
        LEqStComp2.append(EQst_copy)
# This selection returns a list of tuples (staftion_n, earthquake_m)
# Write sequence of lines at the end of the file.ll
json.dump(LEqStComp2,Outputfile)
Outputfile.close()
