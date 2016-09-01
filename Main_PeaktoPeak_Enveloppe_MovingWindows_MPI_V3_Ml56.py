# -*- coding: utf-8 -*-
"""
Created on Tue May 24 16:17:21 2016

@author: claire

Main Peak to Peak Eveloppe
Peak to Peak Maximum of Enveloppe and Integral of enveloppe is calculated for a given station, for all the earthquakes calculated. 
 Main Peak to Peak on moving window 
1. peak to peak ware calculate a save for a selection of events and then save in a file for each station (one ligne per event each column correcpond toa given band pass)
2. Ratio are made compare to a reference station typicaly '5' or '7'
3. Mean of the ratio are calculated and plot
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
MlMin= 5
MlMax = 6
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

with open('/home/rault/PHD/Data/List_EQ_St_Comp_Freq_Ml_%s_%s_distmax_%s.txt'%(MlMin,MlMax,DistMax)) as f:
    LEqStCompFreq =json.load(f)

"""2. creation of a h5py file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""
f = h5py.File('/home/rault/PHD/Results/Parallel_PeaktoPeak_test2_V3_Ml_%s_%s_distmax_%s.hdf5'%(MlMin,MlMax,DistMax), 'w',driver= 'mpio', comm=Comm,libver='latest')
#use libver='latest' will help to be faster !
""" 3. Split the list of Component_frequency and the list of dataSet ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""
def _split_seq(seq, size):

    newseq = []
    splitsize = 1.0/size*len(seq)
    for i in range(size):
            newseq.append(seq[int(round(i*splitsize)):int(round((i+1)*splitsize))])
    return newseq
    
    
dt = dt = h5py.special_dtype(vlen=unicode)
#split the list between the workers
LEqStCompFreq_split = _split_seq(LEqStCompFreq,size)[me]

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
                dStPtP = f[EQname]['St{}'.format(Station)].create_dataset('PtP_{}'.format(Component),(11,), dtype='f')
                dSecondS = f[EQname]['St{}'.format(Station)].create_dataset('SecondS_{}'.format(Component),(11,), dtype='f')
                dSecondP = f[EQname]['St{}'.format(Station)].create_dataset('SecondP_{}'.format(Component),(11,), dtype='f')
                dfreq = f[EQname]['St{}'.format(Station)].create_dataset('Freq_{}'.format(Component),(11,), dtype='f')
                dStPGA = f[EQname]['St{}'.format(Station)].create_dataset('PGA_{}'.format(Component),shape=(1,),dtype='f')
                dStPGV = f[EQname]['St{}'.format(Station)].create_dataset('PGV_{}'.format(Component),shape=(1,),dtype='f')
                dStArias = f[EQname]['St{}'.format(Station)].create_dataset('Arias_{}'.format(Component),shape=(1,),dtype='f')
                dStSatured = f[EQname]['St{}'.format(Station)].create_dataset('Satured_{}'.format(Component),shape=(1,),data=True)
                dStvalid = f[EQname]['St{}'.format(Station)].create_dataset('valid_{}'.format(Component),shape=(1,),data=False)
                dStSNR = f[EQname]['St{}'.format(Station)].create_dataset('SNR_{}'.format(Component),shape=(1,),dtype='f')
                dStart= f[EQname]['St{}'.format(Station)].create_dataset('Start_time_{}'.format(Component),shape=(1,),dtype='f') 
                dSend= f[EQname]['St{}'.format(Station)].create_dataset('End_time_{}'.format(Component),shape=(1,),dtype='f') 
                dStexist = f[EQname]['St{}'.format(Station)].create_dataset('Exist_{}'.format(Component),shape=(1,),data=True)
    if 'St{}'.format(Station) not in f[EQname]: 
            #print 'creation group Station ' , Station , 'EQ ', EQname
            grpSt = f[EQname].create_group('St{}'.format(Station))
            #save the metadata of the Station 
            grpSt.attrs['RDist'] = Rdistance
            grpSt.attrs['Az'] = Az
            grpSt.attrs['BAz'] = BAz
            
            for Component in LComponents:
                if 'PtP_{}'.format(Component) not in f[EQname]['St{}'.format(Station)]:
                    dStPtP = f[EQname]['St{}'.format(Station)].create_dataset('PtP_{}'.format(Component),(11,), dtype='f')
                    dSecondS = f[EQname]['St{}'.format(Station)].create_dataset('SecondS_{}'.format(Component),(11,), dtype='f')
                    dSecondP = f[EQname]['St{}'.format(Station)].create_dataset('SecondP_{}'.format(Component),(11,), dtype='f')
                    dfreq = f[EQname]['St{}'.format(Station)].create_dataset('Freq_{}'.format(Component),(11,), dtype='f')
                    dStPGA = f[EQname]['St{}'.format(Station)].create_dataset('PGA_{}'.format(Component),shape=(1,),dtype='f')
                    dStPGV = f[EQname]['St{}'.format(Station)].create_dataset('PGV_{}'.format(Component),shape=(1,),dtype='f')
                    dStArias = f[EQname]['St{}'.format(Station)].create_dataset('Arias_{}'.format(Component),shape=(1,),dtype='f')
                    dStSatured = f[EQname]['St{}'.format(Station)].create_dataset('Satured_{}'.format(Component),shape=(1,),data=True)
                    dStvalid = f[EQname]['St{}'.format(Station)].create_dataset('valid_{}'.format(Component),shape=(1,),data=False)
                    dStSNR = f[EQname]['St{}'.format(Station)].create_dataset('SNR_{}'.format(Component),shape=(1,),dtype='f')
                    dStart= f[EQname]['St{}'.format(Station)].create_dataset('Start_time_{}'.format(Component),shape=(1,),dtype='f') 
                    dSend= f[EQname]['St{}'.format(Station)].create_dataset('End_time_{}'.format(Component),shape=(1,),dtype='f') 
                    dStexist = f[EQname]['St{}'.format(Station)].create_dataset('Exist_{}'.format(Component),shape=(1,),data=True)
print 'database built'

""" 6. Itaration over EQ, ST, Comp and frequency ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""
L=[]
for worker in range(size):

    if worker==me : 
        for EQstCompFreq in LEqStCompFreq_split:
        #            # 6.0 take case feature EQ, station component and frequency
            EQname, Station, Year,jJul, Hour, Secondp, Seconds,Ml,Depth, Rdistance, Lat,Long, Az,BAz, Component, number, freqband = EQstCompFreq
                 
                
            if Year<>f[EQname].attrs['Year'] or  Hour<>f[EQname].attrs['Hour'] or  jJul<>f[EQname].attrs['JJul']: 
                print 'ERROR PROBLEM OF EQ FEATURES' , 'Je suis le worker {} et je traite le seisme {}, et la composante {} pour la frequence max {}'.format(me,EQname,Component, freqband[1])
                
                
            ## 6.1 read file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if not  os.path.exists('/home/burtin/DATA/LinTianShan/Seismic_Data/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHN.SAC'%(Year,jJul,Station,Year,jJul,str(Hour))) : 
               print 'file not existing Year :',  Year, 'Julian day : ',jJul, 'Hour : ', Hour
               f[EQname]['St{}'.format(Station)]['Exist_{}'.format(Component)][0]=False
               continue
           
            st=Read_event(Year,jJul,Hour,Secondp, Station, False)
            #if the begining of the signal is too close to 00min then the wave peaking won't be accurate thus it's necessary to merge the signal whith the previous one! 
            if Secondp<100:
                if int(Hour)-1>= 0:
                    st = st + Read_event(Year,jJul,str(int(Hour)-1),Secondp, Station, False)
                else : 
                    if int(jJul)-1>1:
                        st = st + Read_event(Year,int(jJul-1),str(23),Secondp, Station, False)
                    else : 
                        st = st + Read_event(str(int(Year)-1),str(365),str(23),Secondp, Station, False)
            st.merge(method= 1)
            
            if Component in ['R', 'T']:
                st_copy = st.copy()
                strot_copy = Rotation(st_copy ,BAz, False)
                tr1 = strot_copy.select(component=Component)[0]

            elif Component=='H':
                st_copy = st.copy()
                strot_copy= strot.copy()
                tr1, trZabs = RealafterrottrH(strot_copy,False)
            else :
                st_copy = st.copy()
                tr1 = st_copy.select(component=Component)[0]

            #6.1 Satureation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Satured = Saturation(tr1,False)
            
            #6.2 Instrument correction~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            stcorrect = Stream_Correction(st, Station, False)
            
            #6.3 Rotation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            strot = Rotation(stcorrect,BAz, False)
            st = stcorrect.copy()
#
            print 'Je suis le worker {} et je traite le seisme {}, et la composante {} pour la frequence max {}'.format(me,EQname,Component, freqband[1])
#           print 'There is still 
            #6.4 select the trace of the component~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if Component in ['R', 'T']:
                st_copy = st.copy()
                strot_copy= strot.copy()
                tr = strot_copy.select(component=Component)[0]
                trZ = st_copy.select(component="Z")[0]
            elif Component=='H':
                st_copy = st.copy()
                strot_copy= strot.copy()
                tr, trZabs = RealafterrottrH(strot_copy,False)
                trZ = st_copy.select(component="Z")[0]
            else :
                st_copy = st.copy()
                tr = st_copy.select(component=Component)[0]
                trZ = st_copy.select(component="Z")[0]

            #6.5Define Python instance nameEQ,magnitude,depth,Rdistance,Lat,Long,Az,time,Station, frequency band, stream, streamrot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            EQ = PtPearthquake(EQname,f[EQname].attrs['Ml'],f[EQname].attrs['Depth'],f[EQname]['St{}'.format(Station)].attrs['RDist'],f[EQname].attrs['Lat'],f[EQname].attrs['Long'],f[EQname]['St{}'.format(Station)].attrs['Az'],[Year,jJul,Hour,str(int(Secondp)),str(int(Seconds))],Station,freqband,tr,Component)         
          
            #6.6 P and S arrivals calculated on the Z componenent~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         
            EQ.second_P, EQ.second_S = WavePicking2(trZ,5,int(EQ.time[3]),int(EQ.time[4]),False)
          
            #6.7 Trim the signal between ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            EQ.trTrim, starttime,endtime =  Trimtr(EQ.tr,EQ.second_P,1.5, 0.65,7,15,False,False)

            #6.8SNR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            EQ.valid,EQ.SNR= SNRstd(EQ.tr, EQ.trTrim, LimSNR, EQ.second_P, False)
           
            #6.9saturation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            EQ.Satured = Satured
           
#            #6.9bisplot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#            fig = plt.figure()
#            if EQ.valid ==True:
#                val = "T"
#            else : 
#                val ="F"
#            if Station=='1' and Component = N and number =10:  
#                plt.figure()
#                plt.plot(EQ.trTrim.data)
#                plt.title (' valid'+val+'Start '+str(EQ.second_P)+ 'EQ '+ EQname  +'component ' + Component + 'st '+Station )
#                plt.savefig('/home/rault/PHD/Results/EQ_'+ EQ.nameEQ  +'_component_' + EQ.Component + '_st_'+EQ.Station+'.png')
#                print 'Snr calcul',EQ.SNR, val, 'valid'
            #6.10 PGA, PGV, Arias, Peak to Peak~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if EQ.valid==True:
                #PGAG
                EQ.FPGA()
                #PGV
                EQ.FPGV()
                #Arias
                EQ.FArias()
                #2.6 .4 Peak to Peak and filter
                EQ.PtP_envelope()
             
##                2.7 Dictionary containing parameters definition of the EQ considered at the given station for a given frequency of filtering
    ##                dataPickle = qstatcPickle.dumps(EQ)
            Dict = EQ.__dict__
            
            f[EQ.nameEQ]['St{}'.format(EQ.Station)]['PtP_{}'.format(EQ.Component)][number] = Dict['PeaktoPeak']
            f[EQ.nameEQ]['St{}'.format(EQ.Station)]['SecondS_{}'.format(EQ.Component)][number] = Dict['second_S']
            f[EQ.nameEQ]['St{}'.format(EQ.Station)]['SecondP_{}'.format(EQ.Component)][number] = Dict['second_P']
            f[EQ.nameEQ]['St{}'.format(EQ.Station)]['Freq_{}'.format(EQ.Component)][number] = Dict['frband'][0]
            if number==10 : 
                print 'NUMBER',number,'satured', Dict['Satured'][0], 'valid',Dict['valid'],Dict['PGA']
                f[EQ.nameEQ]['St{}'.format(EQ.Station)]['valid_{}'.format(EQ.Component)][0]=Dict['valid']
                f[EQ.nameEQ]['St{}'.format(EQ.Station)]['SNR_{}'.format(EQ.Component)][0] = Dict['SNR']
                f[EQ.nameEQ]['St{}'.format(EQ.Station)]['Satured_{}'.format(EQ.Component)][0] = Dict['Satured'][0]
                f[EQ.nameEQ]['St{}'.format(EQ.Station)]['PGA_{}'.format(EQ.Component)][0] = Dict['PGA']
                f[EQ.nameEQ]['St{}'.format(EQ.Station)]['PGV_{}'.format(EQ.Component)][0] = Dict['PGV']
                f[EQ.nameEQ]['St{}'.format(EQ.Station)]['Arias_{}'.format(EQ.Component)][0]= Dict['Arias']
                #f[EQ.nameEQ]['St{}'.format(EQ.Station)]['EQ_Trim_{}'.format(EQ.Component)][0]= Dict['trTrim']
                f[EQ.nameEQ]['St{}'.format(EQ.Station)]['Start_time_{}'.format(EQ.Component)][0]=calendar.timegm(starttime.utctimetuple())
                f[EQ.nameEQ]['St{}'.format(EQ.Station)]['End_time_{}'.format(EQ.Component)][0]=calendar.timegm(endtime.utctimetuple())
                
         
print 'fclosing'          
# Close h5file
f.close()
# Finalize
MPI.Finalize()
# End time
timeend = MPI.Wtime()
print ('---------------------------------------------/n','time start : ',timestart,'/n','time end : ',timeend,'\n---------------------------------------------/n')
