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
if len(sys.argv) > 1:
  ifile = int(sys.argv[1])

filename = filenames[ifile]
T = 100

FREQS=[1.0*i for i in range(0,21)]
  
splitted = filename.split('_')

Caflux = float(splitted[3])
Lflux = float(splitted[4])
Gluflux = float(splitted[5])
AChflux = float(splitted[6])

isGluR = [i for i in range(0,len(splitted)) if splitted[i].find('GluR') > -1]
isk = [i for i in range(0,len(splitted)) if splitted[i][0]=='k']
 
blocked = '_'.join(splitted[min(isGluR):isk[0]])

print('blocked = '+blocked)

altered = '567,568,569,570,571,572,573,574,575,576x0,0,0,0,0,0,0,0,0,0'

for ifreq in range(0,len(FREQS)):
  freq = max(0.001,FREQS[ifreq])
  Nstim = int(T*freq)
  extfilename = 'nrn_'+str(Nstim)+'_'+str(freq)+'_'+str(Caflux)+'_'+str(Lflux)+'_'+str(Gluflux)+'_'+str(AChflux)+'_'+blocked+'_control'

  if not exists(extfilename+'.mat'):
    print('python3 model_nrn_altered_noU_extfilename_smallconcs_withSpontGluRPhos.py 27000000 1e-6 24040000 '+str(Nstim)+' '+str(freq)+' 3.0 '+str(Caflux)+' '+str(Lflux)+' '+str(Gluflux)+' '+str(AChflux)+' 1 1000 None '+blocked.split('x')[0]+' '+blocked.split('x')[1]+' '+altered.split('x')[0]+' '+altered.split('x')[1]+' '+extfilename)
    os.system('python3 model_nrn_altered_noU_extfilename_smallconcs_withSpontGluRPhos.py 27000000 1e-6 24040000 '+str(Nstim)+' '+str(freq)+' 3.0 '+str(Caflux)+' '+str(Lflux)+' '+str(Gluflux)+' '+str(AChflux)+' 1 1000 None '+blocked.split('x')[0]+' '+blocked.split('x')[1]+' '+altered.split('x')[0]+' '+altered.split('x')[1]+' '+extfilename)
  else:
    print(extfilename+'.mat exists')

