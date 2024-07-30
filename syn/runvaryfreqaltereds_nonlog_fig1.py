import scipy.io
from pylab import *
from os.path import exists
import os

filenames = [
  'nrn_20_0.2_100.0_5.0_10.0_10.0_GluR1,GluR1_memb,GluR2,GluR2_memb,Ca,Cax0.5,0.5,1.5,1.5,1.0,1.0_k1x1.0.mat',
  'nrn_20_0.2_100.0_0.2_10.0_10.0_GluR1,GluR1_memb,GluR2,GluR2_memb,PP1,Cax2.0,2.0,0.0,0.0,0.5,1.0_k1x1.0.mat',
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
ksblockeds = [[179,224,243,288],                 #Phosphorylation of S845 made impossible by blocking PKA binding to GluR1
              [167],                             #Phosphorylation of I1 made impossible by blocking PKA binding to I1
              [167,179,224,243,288],             #Phosphorylation of S845 made impossible by blocking PKA binding to GluR1 and that of I1 by blocking PKA binding to I1
              [434,437,440,443,449,452,455,458], #Phosphorylation of S880 made impossible by blocking PKC binding to GluR2
              [191,194,197,200,212,215,218,221,255,258,261,264,276,279,282,285,434,437,440,443,449,452,455,458], #Phosphorylation of S831 and S880 by PKC made impossible
              [182,185,188,203,206,209,246,249,252,267,270,273], #Phosphorylation of S831 by CK made impossible by blocking CK binding to GluR1
              [191,194,197,200,212,215,218,221,255,258,261,264,276,279,282,285], #Phosphorylation of S831 by PKC made impossible by blocking PKC binding to GluR1
              [182,185,188,203,206,209,246,249,252,267,270,273,191,194,197,200,212,215,218,221,255,258,261,264,276,279,282,285]] #Phosphorylation of S831 made impossible by blocking CK and PKC binding to GluR1

ksnumcoeffs = [[1,0,0], [1,0,0], [1,0,0], [0,0,1], [0,0,1], [0,0,0], [0,0,0], [0,0,0]]

filename = filenames[ifile]
T = 100
Caflux = 100.0
if filename.find(',R') > -1 and filename.find('Leak') == -1:
  Caflux = 150.0
Lflux = 5.0
Gluflux = 10.0
AChflux = 10.0

FREQS=[1.0*i for i in range(0,36)]

if exists(filename):
  A=scipy.io.loadmat(filename)
  
  splitted = filename.split('_')
  isGluR = [i for i in range(0,len(splitted)) if splitted[i].find('GluR') > -1]
  isk = [i for i in range(0,len(splitted)) if splitted[i][0]=='k']
 
  blocked = '_'.join(splitted[min(isGluR):isk[0]])

  print('blocked = '+blocked)
  logcoeffs_all = []
  B = scipy.io.loadmat('logcoeffs_'+blocked+'.mat')
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
      freq = FREQS[ifreq]
      Nstim = int(T*freq)
      extfilename = 'nrn_'+str(Nstim)+'_'+str(freq)+'_'+str(Caflux)+'_'+str(Lflux)+'_'+str(Gluflux)+'_'+str(AChflux)+'_'+blocked+'_icase'+str(icase)

      if not exists(extfilename+'.mat'):
        print('python3 model_nrn_altered_noU_extfilename_smallconcs_withSpontGluRPhos.py 27000000 1e-6 24040000 '+str(Nstim)+' '+str(freq)+' 3.0 '+str(Caflux)+' '+str(Lflux)+' '+str(Gluflux)+' '+str(AChflux)+' 1 1000 None '+blocked.split('x')[0]+' '+blocked.split('x')[1]+' '+altered.split('x')[0]+' '+altered.split('x')[1]+' '+extfilename)
        os.system('python3 model_nrn_altered_noU_extfilename_smallconcs_withSpontGluRPhos.py 27000000 1e-6 24040000 '+str(Nstim)+' '+str(freq)+' 3.0 '+str(Caflux)+' '+str(Lflux)+' '+str(Gluflux)+' '+str(AChflux)+' 1 1000 None '+blocked.split('x')[0]+' '+blocked.split('x')[1]+' '+altered.split('x')[0]+' '+altered.split('x')[1]+' '+extfilename)

else:
  print(filename+' does not exist')
