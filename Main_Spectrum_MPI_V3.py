# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 17:35:24 2016

@author: claire
"""

# -*- coding: utf-8 -*-
"""
Created on Tue May 24 16:17:21 2016

@author: claire

Main EQ Spectrum
Smooth spectrum is calculated for all EQ, station component couple. 

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
from Instance_EQSpectrum import *
import cPickle

#secteion of EQ : 
#Magnitude,min and max, and Distance min and max
MlMin= 3
MlMax = 3.1
DistMin = 0
DistMax = 50
jJultoplot = 206
refStation = '5'
LStations = ['1','2','3','4','5','6','7']
# 0. limit SNR above signal is valid
LimSNR = 2

# 1. Selection (in a separate script) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
File = '/home/rault/PHD/Data/Event_stations_Caracteristics_2015_2016Full2.txt' 

# This selection returns a list of tuples (station_n, earthquake_m) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A_St_EQ = Array_St_EQ_MlDist(File, MlMin, MlMax, DistMin, DistMax)
np.savetxt('/home/rault/PHD/Data/A_St_EQ_ML%s_%s_Test.txt'%(MlMax,MlMax),A_St_EQ)

print 'number of  earthquakes selected with Ml %s max distence is %s'%(MlMax,MlMax), A_St_EQ.shape

#2. creation of a h5py file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f = h5py.File('/home/rault/PHD/Results/Parallel_Spectrum_Ml_%s_%s_distmax_%s.hdf5'%(MlMin,MlMax,DistMax),'w',driver= 'mpio', comm=Comm,libver='latest')

# 3. Create all the groups and datasets~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ldata=[]
#create the empty list of Eq and stations where Eq are records (keep the same order to avoid seismes where there is no records for a considered station)
for EqSta in A_St_EQ:
    #Earthquakes metadata (one line correspond to one EQ and one station)
    Station, Year,jJul, Hour, Secondp,Seconds,Ml,Depth, Rdistance, Lat,Long, Az = ConvertDatestr(EqSta)
    BAz = EqSta[9]
    EQname = str(Year)+'_'+str(jJul)+ '_'+str(Hour)+'_'+str(int(Secondp/60))   
    #list of EQ with its parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Ldata.append([EQname, Station, Year,jJul, Hour, Secondp,Seconds,Ml,Depth, Rdistance, Lat,Long, Az,BAz])
    
    # Create EQ/ST group h5py ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
        for component in ['E','N','Z','R','T', 'H']:
            dStSatured = f[EQname]['St{}'.format(Station)].create_dataset('Satured_{}'.format(component),shape=(1,),data=True)
            dStvalid = f[EQname]['St{}'.format(Station)].create_dataset('valid_{}'.format(component),shape=(1,),data=False)
            dStSNR = f[EQname]['St{}'.format(Station)].create_dataset('SNR_{}'.format(component),shape=(1,),dtype='f')
            #create an empty dataset
            dStSPectrum = f[EQname]['St{}'.format(Station)].create_dataset('Spectrum_{}'.format(component),shape=(50299,1),dtype='f')
    if 'St{}'.format(Station) not in f[EQname]: 
        grpSt = f[EQname].create_group('St{}'.format(Station))
        #save the metadata of the Station 
        grpSt.attrs['RDist'] = Rdistance
        grpSt.attrs['Az'] = Az
        grpSt.attrs['BAz'] = BAz
        for component in ['E','N','Z','R','T', 'H']:
            dStSatured = f[EQname]['St{}'.format(Station)].create_dataset('Satured_{}'.format(component),shape=(1,),data=True)
            dStvalid = f[EQname]['St{}'.format(Station)].create_dataset('valid_{}'.format(component),shape=(1,),data=False)
            dStSNR = f[EQname]['St{}'.format(Station)].create_dataset('SNR_{}'.format(component),shape=(1,),dtype='f')
            #create an empty dataset
            dStSPectrum = f[EQname]['St{}'.format(Station)].create_dataset('Spectrum_{}'.format(component),shape=(50299,1),dtype='f')


# 3. Split the list of Component_frequency and the list of dataSet ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def _split_seq(seq, size):

    newseq = []
    splitsize = 1.0/size*len(seq)
    for i in range(size):
            newseq.append(seq[int(round(i*splitsize)):int(round((i+1)*splitsize))])
    return newseq

LComponents = ['E','N','Z','R','T', 'H']
#list of couple frq component 
LEQStcompotuples = []
for j in xrange(len(Ldata)):
    for i in xrange(len(LComponents)):
        LEQStcompotuples.append([Ldata[j],LComponents[i]])
        
LEqStCompoTuplesplit = _split_seq(LEQStcompotuples,size)[me]


# 4. Iterate over the lists~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ~~~~~~

for worker in xrange(size):
    if worker == me :
        for StaEqComp in LEqStCompoTuplesplit: 
  
            # Get event characteristics
            [EQname, Station, Year,jJul, Hour, Secondp,Seconds,Ml,Depth, Rdistance, Lat,Long, Az,BAz],Component = StaEqComp
  
            ## 2.2 read file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if not  os.path.exists('/home/burtin/DATA/LinTianShan/Seismic_Data/20%s/R%s.02/GSW0%s.%s.%s.%s.00.00.BHN.SAC'%(Year,jJul,Station,Year,jJul,str(Hour))) : 
               #print 'file not existing Year :',  Year, 'Julian day : ',jJul, 'Hour : ', Hour
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
            # Instrument correction~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            stcorrect = Stream_Correction(st, Station, False)
            st_copy = stcorrect.copy()
            #select the trace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if Component in ['R', 'T']:
                # Rotation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                strot = Rotation(stcorrect,BAz, False)
                strot_copy = strot.copy()
                tr = strot_copy.select(component=Component)[0]
                trZ = st_copy.select(component="Z")[0]
                trE = st_copy.select(component='E')[0]
            elif Component=='H':
                strot = Rotation(stcorrect,BAz, False)
                strot_copy = strot.copy()
                tr, trZabs = RealafterrottrH(strot_copy,plot)
                trZ = st_copy.select(component="Z")[0]
                trE = st_copy.select(component='E')[0]
            else :
                st_copy = stcorrect.copy()
                tr = st_copy.select(component=Component)[0]
                trZ = st_copy.select(component='Z')[0]
                trE = st_copy.select(component='E')[0]
            # Define Python instance nameEQ,magnitude,depth,Rdistance,Lat,Long,Az,time,Station, frequency band, stream, streamrot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
            EQ = Spectrumearthquake(EQname,f[EQname].attrs['Ml'],f[EQname].attrs['Depth'],f[EQname]['St{}'.format(Station)].attrs['RDist'],f[EQname].attrs['Lat'],f[EQname].attrs['Long'],f[EQname]['St{}'.format(Station)].attrs['Az'],[Year,jJul,Hour,str(int(Secondp)),str(int(Seconds))],Station,tr,Component)         
            
            # 2.6 P and S arrivals P arrival calculated on the Z and S arrival on the E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            EQ.second_P, second_S = WavePicking2(trZ,5,int(EQ.time[3]),int(EQ.time[4]),False)
            second_P, EQ.second_S = WavePicking2(trE,5,int(EQ.time[3]),int(EQ.time[4]),False)
#            
#            f[EQname]['St{}'.format(Station)]['SecondarrivedP'] = EQ.second_P
#            f[EQname]['St{}'.format(Station)]['SecondarrivedS'] = EQ.second_S
            print 'picking done p start ', EQ.second_P, EQ.second_S
            # 2.6 Trim trace between with a LTA lim = 0.5  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            EQ.trTrim =  Trimtr(EQ.tr,EQ.second_P,1.5, 0.65,7,15,False,False)
            print "trimdone"
            # 2.6.SNR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            EQ.valid,EQ.SNR= SNRstd(EQ.tr, EQ.trTrim, LimSNR, EQ.second_P, False)
            print "snr done",SNR, EQ.valid
            #saturation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            EQ.Satured = Saturation(EQ.trTrim,False)
            print "saturation done"
            if EQ.valid ==True:
                val = "T"
            else : 
                val ="F"
            print "plot"
            plt.figure()
            plt.plot(EQ.trTrim.data)
            plt.title ('valid'+val+'Start '+str(EQ.second_P)+ 'EQ '+ EQname  +'component ' + Component + 'st '+Station )
            plt.savefig('/home/rault/PHD/Results/Spectrum/Spectrum_plotsTraces/EQ'+ EQname  +'component' + Component + 'st'+Station+'.png')
            print 'Snr calcul',EQ.SNR, val, 'valid'
            if EQ.valid==True:
                                            
                #2.6 .0 Spectrum~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                EQ.spectrum = EQ.Spectrum() 
                print "spectrum"
                
                #2.7 Dictionary containing parameters definition of the EQ considered ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            EQ.Save() 
            Dict = EQ.__dict__            
            print 'saving'
            f[EQname]['St{}'.format(Station)]['Spectrum_{}'.format(Component)]['Spectrum_{}'.format(component)]= Dict['spectrum']
            f[EQname]['St{}'.format(Station)]['Spectrum_{}'.format(Component)].attrs['Satured'] = Dict['Satured']
            f[EQname]['St{}'.format(Station)]['Spectrum_{}'.format(Component)].attrs['Valid'] = Dict['valid']
            f[EQname]['St{}'.format(Station)]['Spectrum_{}'.format(Component)].attrs['SNR'] = Dict['SNR']    
            
##2.8 Save in the dataset ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#f[EQname]['St{}'.format(Station)]['Spectrum_{}'.format(Component)]= Dict['spectrum']
#f[EQname]['St{}'.format(Station)]['PtP_{}'.format(Component)].attrs['Satured'] = Dict['Satured']
#f[EQname]['St{}'.format(Station)]['PtP_{}'.format(Component)].attrs['Valid'] = Dict['valid']
#f[EQname]['St{}'.format(Station)]['PtP_{}'.format(Component)].attrs['SNR'] = Dict['SNR']
# # Create EQ/ST group h5py ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    if EQname not in f :   
#        grpEQ = f.create_group(EQname)
#        grpEQ.attrs['Long'] = Long
#        grpEQ.attrs['Lat'] = Lat
#        grpEQ.attrs['Ml'] = Ml
#        grpEQ.attrs['Depth'] = Depth
#        grpEQ.attrs['Year'] = Year 
#        grpEQ.attrs['JJul'] = jJul 
#        grpEQ.attrs['Hour'] = Hour 
#        grpEQ.attrs['Second'] = Second  
#        #create a groupe for station (each line correspond to a EQ station couple)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        grpSt = f[EQname].create_group('St{}'.format(Station))
#        #save the metadata of the Station 
#        grpSt.attrs['RDist'] = Rdistance
#        grpSt.attrs['Az'] = Az
#        grpSt.attrs['BAz'] = BAz
#        #print 'group done ', Station, EQname
#        
#    if 'St{}'.format(Station) not in f[EQname]: 
#            #print 'creation group Station ' , Station , 'EQ ', EQname
#            grpSt = f[EQname].create_group('St{}'.format(Station))
#            #save the metadata of the Station 
#            grpSt.attrs['RDist'] = Rdistance
#            grpSt.attrs['Az'] = Az
#            grpSt.attrs['BAz'] = BAz
#   
#    # Iterate over the workers has to calculate a couple component, frequency band 
#    if 'Spectrum_{}'.format(Component) not in f[EQname]['St{}'.format(Station)]:
#        arr = np.arange(100)
#        dStPtP = f[EQname]['St{}'.format(Station)].create_dataset('Spectrum_{}'.format(Component),data=arr, dtype='i')
#        dStPtP.attrs['Satured'] = 0
#        dStPtP.attrs['SNR'] = 0
#        dStPtP.attrs['Valid'] = 0
            

#print 'fclosing'          
## Close h5file
f.close()

# Finalize
MPI.Finalize()

# End time
timeend = MPI.Wtime()
print ('---------------------------------------------/n','time start : ',timestart,'/n','time end : ',timeend,'\n---------------------------------------------/n')
