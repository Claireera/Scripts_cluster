# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 11:55:39 2016

@author: claire
"""
import json

MlMin= 3
MlMax = 3.1
with open('/home/claire/PHD/Working/Results/List_EQ_St_Comp_Freq_Ml_%s_%s.txt'%(MlMax,MlMax)) as f:
    LEqStCompFreq =json.load(f)
L=[]
i=0
EQ = LEqStCompFreq[462][0]
for line in LEqStCompFreq[462:len(LEqStCompFreq)]:
    print EQ, line[0]
    if line[0]==EQ and line[1]=='1':
        
        L.append(line) 

Outputfile = open('/home/claire/PHD/Working/Data/List_EQ_St_Comp_Freq_Ml_%s_%sTEST.txt'%(MlMax,MlMax), mode = 'w+')
json.dump(L,Outputfile)
Outputfile.close()
#
#with open('/home/claire/PHD/Working/Data/List_EQ_St_Comp_Freq_Ml_%s_%sTEST.txt'%(MlMax,MlMax)) as f:
#    LEqStCompFreq =json.load(f)
#    L7 = []
#    for line in LEqStCompFreq :
#        if line[1]=='3':
#            L7.append(line)