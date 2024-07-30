import scipy.io
from pylab import *
from os.path import exists
import os

filenames = [
  'nrn_20_0.2_150.0_5.0_10.0_10.0_GluR1,GluR1_memb,GluR2,GluR2_memb,Ca,Cax0.5,0.5,1.5,1.5,1.0,1.0_k1x1.0.mat',
  'nrn_20_0.2_150.0_0.2_10.0_10.0_GluR1,GluR1_memb,GluR2,GluR2_memb,PP1,Cax2.0,2.0,0.0,0.0,0.5,1.0_k1x1.0.mat',
  'nrn_20_0.2_150.0_5.0_10.0_10.0_GluR1,GluR1_memb,GluR2,GluR2_memb,R,PP1,Cax2.0,2.0,0.0,0.0,0.0,0.0,1.0_k1x1.0.mat'
]


ifile = 0
icasegiven = -1
if len(sys.argv) > 1:
  ifile = int(sys.argv[1])
if len(sys.argv) > 2:
  icasegiven = int(sys.argv[2])

ksnums = [[567,568,569,570], #S845
          [571,572,573,574], #S831
          [575,576]]         #S880
#ksblockeds = [[373,374], [416,417,418], [416,417,418]] #Only spont. 845 when PKA blocked, only spont. 880 when PKC blocked, only spont. 831 when PKC blocked
ksblockeds = [[181,226,245,290],                 #Phosphorylation of S845 made impossible by blocking phosphorylation of PKAc-bound GluR1                                
              [169],                             #Phosphorylation of I1 made impossible by blocking phosphorylatoin of PKAc-bound I1                                     
              [169,181,226,245,290],             #Phosphorylation of S845 made impossible by blocking PKA-mediated phosphorylation of GluR1 and I1                       
              [436,439,442,445,451,454,457,460], #Phosphorylation of S880 made impossible by blocking phosphorylation of PKC-bound GluR2                                 
              [193,196,199,202,214,217,220,223,257,260,263,266,278,281,284,287,436,439,442,445,451,454,457,460], #Phosphorylation of S831 and S880 by PKC made impossible
              [184,187,190,205,208,211,248,251,254,269,272,275], #Phosphorylation of S831 by CK made impossible
              [193,196,199,202,214,217,220,223,257,260,263,266,278,281,284,287], #Phosphorylation of S831 by PKC made impossible
              [184,187,190,205,208,211,248,251,254,269,272,275,193,196,199,202,214,217,220,223,257,260,263,266,278,281,284,287]] #Phosphorylation of S831 by CK and PKC made impossible


ksnumcoeffs = [[1,0,0], [1,0,0], [1,0,0], [0,0,1], [0,0,1], [0,0,0], [0,0,0], [0,0,0]]

filename = filenames[ifile]
T = 100

FREQS=[1.0*i for i in range(0,21)]

if True: #exists(filename):
  #A=scipy.io.loadmat(filename) #this is not needed, is it?
  
  splitted = filename.split('_')

  Caflux = float(splitted[3])
  Lflux = float(splitted[4])
  Gluflux = float(splitted[5])
  AChflux = float(splitted[6])

  isGluR = [i for i in range(0,len(splitted)) if splitted[i].find('GluR') > -1]
  isk = [i for i in range(0,len(splitted)) if splitted[i][0]=='k']
 
  blocked = '_'.join(splitted[min(isGluR):isk[0]])

  print('blocked = '+blocked)
  logcoeffs_all = []
  print('Loading logcoeffs2_'+blocked+'.mat')
  B = scipy.io.loadmat('logcoeffs2_'+blocked+'.mat')
  logcoeffs_all = B['logcoeffs_all']

  for icase in [0,3]:
    if icasegiven > -1 and icase != icasegiven:
      continue
    if icase < len(logcoeffs_all):
      logcoeffs = logcoeffs_all[icase]
    else:
      logcoeffs = [-20,-20,-20] #set all spontaneous rate coefficients to zero if these are extra runs
    ksblocked = ksblockeds[icase]
    ksnumcoeff = ksnumcoeffs[icase]

    kcoeff = 10**logcoeffs[2]

    altered = ','.join([str(x) for x in ksblocked])+',567,568,569,570,571,572,573,574,575,576'+'x'+\
              ','.join(['0' for x in ksblocked])+','+','.join([','.join([str(kcoeff*ksnumcoeffs[icase][ik]) for i in ksnums[ik]]) for ik in [0,1,2]])

    for ifreq in range(0,len(FREQS)):
      freq = max(0.001,FREQS[ifreq])
      Nstim = int(T*freq)
      extfilename = 'nrn_altered2_'+str(Nstim)+'_'+str(freq)+'_'+str(Caflux)+'_'+str(Lflux)+'_'+str(Gluflux)+'_'+str(AChflux)+'_'+blocked+'_icase'+str(icase)

      if not exists(extfilename+'.mat'):
        print('python3 model_nrn_altered_noU_extfilename_smallconcs_withSpontGluRPhos.py 27000000 1e-6 24040000 '+str(Nstim)+' '+str(freq)+' 3.0 '+str(Caflux)+' '+str(Lflux)+' '+str(Gluflux)+' '+str(AChflux)+' 1 1000 None '+blocked.split('x')[0]+' '+blocked.split('x')[1]+' '+altered.split('x')[0]+' '+altered.split('x')[1]+' '+extfilename)
        os.system('python3 model_nrn_altered_noU_extfilename_smallconcs_withSpontGluRPhos.py 27000000 1e-6 24040000 '+str(Nstim)+' '+str(freq)+' 3.0 '+str(Caflux)+' '+str(Lflux)+' '+str(Gluflux)+' '+str(AChflux)+' 1 1000 None '+blocked.split('x')[0]+' '+blocked.split('x')[1]+' '+altered.split('x')[0]+' '+altered.split('x')[1]+' '+extfilename)
      else:
        print(extfilename+'.mat exists')

else:
  print(filename+' does not exist')
