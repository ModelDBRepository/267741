#cp printbaseline_doiters.py printbaseline_doiters2.py #changed alteration of binding to GluR to that of phosphorylation
#cp printbaselineS831_S845_S880_doiterWithSpontGluRPhos.py printbaseline_doiters.py
import scipy.io
from pylab import *
from os.path import exists
import os

filename = sys.argv[1]
allowExisting = 0
if len(sys.argv) > 2:
  allowExisting = int(sys.argv[2])

ksnums = [[567,568,569,570], #S845
          [571,572,573,574], #S831
          [575,576]]         #S880

ksblockeds = [[181,226,245,290],                 #Phosphorylation of S845 made impossible by blocking phosphorylation of PKAc-bound GluR1
              [169],                             #Phosphorylation of I1 made impossible by blocking phosphorylatoin of PKAc-bound I1
              [169,181,226,245,290],             #Phosphorylation of S845 made impossible by blocking PKA-mediated phosphorylation of GluR1 and I1
              [436,439,442,445,451,454,457,460], #Phosphorylation of S880 made impossible by blocking phosphorylation of PKC-bound GluR2
              [193,196,199,202,214,217,220,223,257,260,263,266,278,281,284,287,436,439,442,445,451,454,457,460]] #Phosphorylation of S831 and S880 by PKC made impossible


ksnumcoeffs = [[1,0,0], [1,0,0], [1,0,0], [0,0,1], [0,0,1]]
itargets = [3,3,3,4,4]


if exists(filename):
  A=scipy.io.loadmat(filename)
  iS831 = [i for i in range(0,len(A['headers'])) if 'S831' in A['headers'][i] ]
  iS845 = [i for i in range(0,len(A['headers'])) if 'S845' in A['headers'][i] ]
  iS880 = [i for i in range(0,len(A['headers'])) if 'S880' in A['headers'][i] ]
  iGluR1memb = [i for i in range(0,len(A['headers'])) if 'GluR1_memb' in A['headers'][i] ]
  iGluR2memb = [i for i in range(0,len(A['headers'])) if 'GluR2_memb' in A['headers'][i] ]
  targetS831 = 1e6*sum([A['DATA'][i,0] for i in iS831])
  targetS845 = 1e6*sum([A['DATA'][i,0] for i in iS845])
  targetS880 = 1e6*sum([A['DATA'][i,0] for i in iS880])
  targetGluR1memb = 1e6*sum([A['DATA'][i,0] for i in iGluR1memb])
  targetGluR2memb = 1e6*sum([A['DATA'][i,0] for i in iGluR2memb])
  print("Baseline S831: "+str(targetS831))
  print("Baseline S845: "+str(targetS845))
  print("Baseline S880: "+str(targetS880))
  print("Baseline GluR1memb: "+str(targetGluR1memb))
  print("Baseline GluR2memb: "+str(targetGluR2memb))
  targets = [targetS831,targetS845,targetS880,targetGluR1memb,targetGluR2memb]
  
  splitted = filename.split('_')
  isGluR = [i for i in range(0,len(splitted)) if splitted[i].find('GluR') > -1]
  isk = [i for i in range(0,len(splitted)) if splitted[i][0]=='k']
 
  blocked = '_'.join(splitted[min(isGluR):isk[0]])

  print('blocked = '+blocked)
  logcoeffs_all = []
  for icase in [0,3]: #cases 1 and 2 and 4 not needed.
    ksblocked = ksblockeds[icase]
    ksnumcoeff = ksnumcoeffs[icase]
    logcoeffs = [-9,0,-4.5]
    myvals_considered = []
    Niter = 20
    for iter in range(0,Niter):
      print(str(logcoeffs))
      kcoeff = 10**logcoeffs[min(iter,2)]
      altered = ','.join([str(x) for x in ksblocked])+',567,568,569,570,571,572,573,574,575,576'+'x'+\
                ','.join(['0' for x in ksblocked])+','+','.join([','.join([str(kcoeff*ksnumcoeffs[icase][ik]) for i in ksnums[ik]]) for ik in [0,1,2]])
      
      extfilename = 'tmp_'+blocked+'_icase'+str(icase)+'_iter'+str(iter)
      if not exists(extfilename+'.mat') or not allowExisting:
        print('python3 model_nrn_altered_noU_extfilename_smallconcs_withSpontGluRPhos.py 25000000 1e-6 24040000 1 1.0 0.0 0.0 0.0 0.0 0.0 1 1000 None '+blocked.split('x')[0]+' '+blocked.split('x')[1]+' '+altered.split('x')[0]+' '+altered.split('x')[1]+' '+extfilename)
        os.system('python3 model_nrn_altered_noU_extfilename_smallconcs_withSpontGluRPhos.py 25000000 1e-6 24040000 1 1.0 0.0 0.0 0.0 0.0 0.0 1 1000 None '+blocked.split('x')[0]+' '+blocked.split('x')[1]+' '+altered.split('x')[0]+' '+altered.split('x')[1]+' '+extfilename)
      A = scipy.io.loadmat(extfilename+'.mat')
      iS831 = [i for i in range(0,len(A['headers'])) if 'S831' in A['headers'][i] ]
      iS845 = [i for i in range(0,len(A['headers'])) if 'S845' in A['headers'][i] ]
      iS880 = [i for i in range(0,len(A['headers'])) if 'S880' in A['headers'][i] ]
      iGluR1memb = [i for i in range(0,len(A['headers'])) if 'GluR1_memb' in A['headers'][i] ]
      iGluR2memb = [i for i in range(0,len(A['headers'])) if 'GluR2_memb' in A['headers'][i] ]
      myS831 = 1e6*sum([A['DATA'][i,0] for i in iS831])
      myS845 = 1e6*sum([A['DATA'][i,0] for i in iS845])
      myS880 = 1e6*sum([A['DATA'][i,0] for i in iS880])
      myGluR1memb = 1e6*sum([A['DATA'][i,0] for i in iGluR1memb])
      myGluR2memb = 1e6*sum([A['DATA'][i,0] for i in iGluR2memb])
      myvals = [myS831, myS845, myS880, myGluR1memb, myGluR2memb]
      myvals_considered.append(myvals[itargets[icase]])
      print("target = "+str(targets[itargets[icase]])+", val = "+str(myvals[itargets[icase]]))
      if iter == 2 and myvals_considered[0] > myvals_considered[1]:
        logcoeffs = [logcoeffs[1], logcoeffs[0], logcoeffs[2]] # this is only done max. 1 time: if it's a descending curve, switch the first two values
        myvals_considered = [myvals_considered[1], myvals_considered[0], myvals_considered[2]]
      if iter > 1 and myvals[itargets[icase]] > targets[itargets[icase]]:
        logcoeffs = [logcoeffs[0],logcoeffs[2], 0.5*(logcoeffs[0]+logcoeffs[2])]
      elif iter > 1:
        logcoeffs = [logcoeffs[2],logcoeffs[1], 0.5*(logcoeffs[2]+logcoeffs[1])]
    logcoeffs_all.append(logcoeffs[:])
  scipy.io.savemat('logcoeffs2_'+blocked+'.mat', {'logcoeffs_all': logcoeffs_all, 'blocked': blocked}) 
  print('logcoeffs_all = '+str(logcoeffs_all))
else:
  print(filename+' does not exist')
