#cp fitsheu_check_varEGTA_feweralt.py fitsheu_check_EGTA2.py
from pylab import *
import os
from os.path import exists
import random
import scipy.io
import emoo
import pickle

#These give for some reason different [Ng] than in the paper, wonder why...
#daltonInug = 1.66054e-18
#conc_Ng = 40.0/(7600*daltonInug)*6.022e-23*1e6 # 40ug/ml -> 40.0/(7600*daltonInug) molecules/ml -> 40.0/(7600*daltonInug)*6.022e-23 mol/ml -> 40.0/(7600*daltonInug)*6.022e-23*1e6 mM
conc_Ng = 0.0053 #mM, from Sheu et al. 1995 p. 337.
conc_PKC = 0.00000857 #mM, from Sheu et al.: 0.66 ug/ml
conc_Ca = 0.1    #mM, when no EGTA used
conc_Ca_EGTA = 0 #mM, when 2mm EGTA used

DATA_Fig6a = [[0.0, 100.0],
  [1.4565826330532197e-3, 83.90243902439023],
  [2.4649859943977592e-3, 42.439024390243915],
  [4.817927170868348e-3, 17.5609756097561],
  [9.635854341736696e-3, 11.21951219512195],
  [19.15966386554622e-3, 11.21951219512195],
  [38.31932773109244e-3, 8.78048780487805]]

DATA_Fig6b = [[0.0, 100.0],
  [1.461684316431695e-3, 93.59190991525531],
  [2.6819043688005486e-3, 73.48665626168943],
  [5.031537399873592e-3, 52.91188876133475],
  [9.761244469668632e-3, 38.70264552994439],
  [18.892514865788698e-3, 19.876946096198765],
  [38.4351258271312e-3, 15.087646884311312]]

myseed = 1
N_samples = 1000
N_check = 1
if len(sys.argv) > 1:
  myseed = int(float(sys.argv[1]))
if len(sys.argv) > 2:
  N_samples = int(float(sys.argv[2]))
if len(sys.argv) > 4:
  N_check = int(float(sys.argv[4]))

def phosphs(ks,EGTAcoeff=0.0,verbose=False,rankID=0):
  NgpCaM_bind_coeff = 1.0
  NgpCaM_unbind_coeff = 1.0

  NgPKC_coeff = ks[0]            #Multiply the backward rate of Ng + PKCt <-> NgPKCt by this coeff (and change forward and catalysis rates to keep the original Vmax and Km)
  NgCaM_PKC_bind_coeff = 1.0     #Multiply the forward rate of NgCaM* + PKCt <-> NgPKCt by this coeff (compared to rate of Ng + PKCt <-> NgPKCt). This will be done after the multiplying of ks[0] (Vmax and Km do not hold for NgCaM*)
  NgCaM_PKC_unbind_coeff = 1.0   #Multiply the backward rate of NgCaM* + PKCt <-> NgPKCt by this coeff (compared to rate of Ng + PKCt <-> NgPKCt). This will be done after the multiplying of ks[0] (Vmax and Km do not hold for NgCaM*)
  NgCaM_PKC_cat_coeff = ks[1]    #Multiply the catalysis rate of NgCaM* + PKCt <-> NgPKCt by this coeff (compared to rate of Ng + PKCt <-> NgPKCt). This will be done after the multiplying of ks[0] (Vmax and Km do not hold for NgCaM*)

  Ng_dephosph = 0.0

  Ngps_tc_all = []
  Ngps_all = []
  Ngps_10min_all = []
  lastConcs_all = []
  headers = []

  for iCaAmount in [0,1]:
    Ngps = []
    Ngps_10min = []
    if iCaAmount == 0:
      CaCoeff = EGTAcoeff
      concs_CaM = [x[0] for x in  DATA_Fig6a]
    else: 
      CaCoeff = 1.0
      concs_CaM = [x[0] for x in  DATA_Fig6b]
    randomString = ''.join(random.choice('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789') for _ in range(0,10))+'_'+str(rankID)
    filenames = []
    for iconc_CaM in range(0,len(concs_CaM)):
      conc_CaM = concs_CaM[iconc_CaM]
      filename = 'fit'+randomString+'_'+str(iconc_CaM)+'.mat'

      filename_IC = 'nrnICs/nrnICs_short_tstop0.0001_Ca-CaM-Ng-PKCt-PKCp-PKCx'+str(1e6*conc_Ca*CaCoeff)+'-'+str(1e6*conc_CaM)+'-'+str(1e6*conc_Ng)+'-'+str(0.5e6*conc_PKC)+'-'+str(0.5e6*conc_PKC)+'-0.0.mat'

      #Some manual work: str() worked differently in python2 vs python3. Here corrected one by one
      for inum in range(0,9):
        filename_IC = filename_IC.replace(str(inum)+'999999999999',str(inum+1))
      filename_IC = filename_IC.replace('5826330532198','58263305')
      filename_IC = filename_IC.replace('9859943977594','9859944')
      filename_IC = filename_IC.replace('927170868348','92717087')
      filename_IC = filename_IC.replace('854341736696','85434174')
      filename_IC = filename_IC.replace('66386554622','6638655')
      filename_IC = filename_IC.replace('32773109244','3277311')
      filename_IC = filename_IC.replace('6843164316952','68431643')
      filename_IC = filename_IC.replace('9043688005486','9043688')
      filename_IC = filename_IC.replace('537399873591','53739987')
      filename_IC = filename_IC.replace('244469668633','24446967')
      filename_IC = filename_IC.replace('514865788697','5148658')
      filename_IC = filename_IC.replace('125827131196','1258271')

      #Load the ICs from filename_IC. Remember that these will be in mM.
      ICMAT = scipy.io.loadmat(filename_IC)
      mystrIC_species = ''
      mystrIC_concs = ''
      for ispec2 in range(0,len(ICMAT['headers'])):
        header = ICMAT['headers'][ispec2]
        if header.find(' ') > -1:
          header = header[0:header.find(' ')]
        mystrIC_species = mystrIC_species + header + ','
        mystrIC_concs = mystrIC_concs + str(ICMAT['DATA'][ispec2][-1]) + ',' #These are printed in mM, like the data in ICMAT
      mystrIC_species = mystrIC_species[0:-1]
      mystrIC_concs = mystrIC_concs[0:-1]

      catperbackw = 887.2520420070011 # (kcat/kbackw = 1.7745040840140023 (1/ms) /0.002 (1/ms)
      NgPKC_coeff_forw = (NgPKC_coeff + catperbackw)/(1.0 + catperbackw)
      if verbose:
        print(('python3 model_nrn_CaM_Ng_PKC_only_altered_extfilename_absconcs.py 3600000 1e-7 '+mystrIC_species+' '+mystrIC_concs+' '+\
          '20,24,28,32,21,25,29,33,36,51,37,52,39,41,43,45,54,56,58,60,40,42,44,46,55,57,59,61,47,48,49,50,62,63,64,65,66,67,68,69,70 '+\
           ','.join([str(NgpCaM_bind_coeff) for i in [0,1,2,3]])+','+','.join([str(NgpCaM_unbind_coeff) for i in [0,1,2,3]])+','+\
           ','.join([str(NgPKC_coeff_forw) for i in [0,1]])+','+','.join([str(NgPKC_coeff) for i in [0,1]])+','+\
           ','.join([str(NgCaM_PKC_bind_coeff*NgPKC_coeff_forw) for i in range(0,8)])+','+','.join([str(NgCaM_PKC_unbind_coeff*NgPKC_coeff) for i in range(0,8)])+','+','.join([str(NgCaM_PKC_cat_coeff) for i in range(0,8)])+','+\
           ','.join([str(Ng_dephosph) for i in range(0,5)])+' '+filename))
        os.system('python3 model_nrn_CaM_Ng_PKC_only_altered_extfilename_absconcs.py 3600000 1e-7 '+mystrIC_species+' '+mystrIC_concs+' '+\
          '20,24,28,32,21,25,29,33,36,51,37,52,39,41,43,45,54,56,58,60,40,42,44,46,55,57,59,61,47,48,49,50,62,63,64,65,66,67,68,69,70 '+\
           ','.join([str(NgpCaM_bind_coeff) for i in [0,1,2,3]])+','+','.join([str(NgpCaM_unbind_coeff) for i in [0,1,2,3]])+','+\
           ','.join([str(NgPKC_coeff_forw) for i in [0,1]])+','+','.join([str(NgPKC_coeff) for i in [0,1]])+','+\
           ','.join([str(NgCaM_PKC_bind_coeff*NgPKC_coeff_forw) for i in range(0,8)])+','+','.join([str(NgCaM_PKC_unbind_coeff*NgPKC_coeff) for i in range(0,8)])+','+','.join([str(NgCaM_PKC_cat_coeff) for i in range(0,8)])+','+\
           ','.join([str(Ng_dephosph) for i in range(0,5)])+' '+filename)
      else:
        os.system('python3 model_nrn_CaM_Ng_PKC_only_altered_extfilename_absconcs.py 3600000 1e-7 '+mystrIC_species+' '+mystrIC_concs+' '+\
          '20,24,28,32,21,25,29,33,36,51,37,52,39,41,43,45,54,56,58,60,40,42,44,46,55,57,59,61,47,48,49,50,62,63,64,65,66,67,68,69,70 '+\
           ','.join([str(NgpCaM_bind_coeff) for i in [0,1,2,3]])+','+','.join([str(NgpCaM_unbind_coeff) for i in [0,1,2,3]])+','+\
           ','.join([str(NgPKC_coeff_forw) for i in [0,1]])+','+','.join([str(NgPKC_coeff) for i in [0,1]])+','+\
           ','.join([str(NgCaM_PKC_bind_coeff*NgPKC_coeff_forw) for i in range(0,8)])+','+','.join([str(NgCaM_PKC_unbind_coeff*NgPKC_coeff) for i in range(0,8)])+','+','.join([str(NgCaM_PKC_cat_coeff) for i in range(0,8)])+','+\
           ','.join([str(Ng_dephosph) for i in range(0,5)])+' '+filename+' > /dev/null')
        #'20-24-28-32-' #Ngp binding to CaM*;     multiply by NgpCaM_bind_coeff
        #'21-25-29-33-' #Ngp unbinding from CaM*; multiply by NgpCaM_unbind_coeff; (this is alreaday tentatively 10 times faster than Ng unbinding from CaM*)
        #'36-51-' #Ng bind to PKC;     multiply by NgPKC_coeff_forw (ensures the right Km)
        #'37-52-' #Ng unbind from PKC; multiply by NgPKC_coeff
        #'38-53-' #NgPKC -> Ngp;       multiply by 1.0 (these rates are determined by Vmax and [PKCt], don't change them here)
        #'39-41-43-45-54-56-58-60-' #NgCaM bind to PKC;     multiply by NgCaM_PKC_bind_coeff*NgPKC_coeff_forw (correct Km if NgCaM_PKC_bind_coeff==1.0)
        #'40-42-44-46-55-57-59-61-' #NgCaM unbind from PKC; multiply by NgCaM_PKC_unbind_coeff*NgPKC_coeff    (correct Km if NgCaM_PKC_unbind_coeff==1.0)
        #'47-48-49-50-62-63-64-65-' #NgCaMPKC catalysis;    multiply by NgCaM_PKC_cat_coeff                   (correct Km if NgCaM_PKC_cat_coeff==1.0)
        #'66'                       # Ng_dephosph
      if not exists(filename):
        print(filename+' does not exist')
        Ngps.append(nan)
        Ngps_10min.append(nan)
        break
      MAT = scipy.io.loadmat(filename)
      Ngp = sum(array([MAT['DATA'][j,:] for j in range(0,len(MAT['headers'])) if 'Ngp' in MAT['headers'][j]]),axis=0)
      Ngps.append(Ngp[-1]/conc_Ng)
      Ngps_tc_all.append(Ngp[:]/conc_Ng)
      NgCaM = sum(array([MAT['DATA'][j,:] for j in range(0,len(MAT['headers'])) if 'Ng' in MAT['headers'][j] and 'CaM' in MAT['headers'][j]]),axis=0)
      filenames.append(filename)
      i10min = argmin([abs(MAT['DATA'][0][j]-600000) for j in range(0,len(MAT['DATA'][0]))])
      Ngps_10min.append(Ngp[i10min]/conc_Ng)
      lastConcs_all.append(MAT['DATA'][:,-1])
      headers = MAT['headers']

    for filename in filenames:
      os.system('rm '+filename)
    while len(Ngps) < len(concs_CaM):
      Ngps.append(nan)
      Ngps_10min.append(nan)
    Ngps_all.append(Ngps[:])
    Ngps_10min_all.append(Ngps_10min[:])
    if isnan(Ngps[-1]):
      if iCaAmount == 0: # If any simulation produced nan, all following simulation results will be nan as well, and the function evaluation is terminated
        Ngps_all.append(Ngps[:])
        Ngps_10min_all.append(Ngps_10min[:])
      break

  return [Ngps_all, Ngps_10min_all, Ngps_tc_all, lastConcs_all, headers]

OBJECTIVES = ['f0','f1','f2','f3']
VARIABLES = [['k0',-7.,2.], #Note: these are log10 values, i.e., each variable is allowed to vary from 0.0000001-fold to 100-fold
             ['k1',-7.,2.]]
MAXERR = 1e8 # Maximum error for missing data

def func_to_optimize(parameters,EGTAcoeff=0.0,saveFig="",rankID=0,verbose=False):
  mydict = {}
  phosphs_nonnorm_both = phosphs([10**parameters[VARIABLES[i][0]] for i in range(0,len(VARIABLES))], EGTAcoeff, verbose, rankID=rankID)
  phosphs_nonnorm = phosphs_nonnorm_both[0]
  phosphs_nonnorm_10min = phosphs_nonnorm_both[1]
  phosphs_nonnorm_tc = phosphs_nonnorm_both[2]
  lastConcs_all = phosphs_nonnorm_both[3]
  headers = phosphs_nonnorm_both[4]
  

  mydict['f0'] = sum([abs(phosphs_nonnorm[0][i]/phosphs_nonnorm[0][0]*100 - DATA_Fig6a[i][1]) for i in range(0,len(DATA_Fig6a))])
  mydict['f1'] = sum([abs(phosphs_nonnorm[1][i]/phosphs_nonnorm[1][0]*100 - DATA_Fig6b[i][1]) for i in range(0,len(DATA_Fig6b))])
  mydict['f2'] = abs(phosphs_nonnorm_10min[0][0]/phosphs_nonnorm[0][0]*100 - 50)
  for i in range(0,2):
    if isnan(mydict['f'+str(i)]):
      mydict['f'+str(i)] = MAXERR

  if len(saveFig) > 0:
    f,axarr = subplots(2,1)
    for iax in range(0,2):
      axarr[iax].set_position([0.1+0.37*iax,0.8,0.25,0.15])
      axarr[iax].set_ylim([0,100])
      f.text(0.03+0.37*iax,0.94,chr(ord('D')+iax),fontsize=16)
      for tick in axarr[iax].xaxis.get_major_ticks() + axarr[iax].yaxis.get_major_ticks():
        tick.label.set_fontsize(5)
      for axis in ['top','bottom','left','right']:
        axarr[iax].spines[axis].set_linewidth(0.5)
        
    axarr[0].plot([DATA_Fig6a[i][0] for i in range(0,len(DATA_Fig6a))],[DATA_Fig6a[i][1] for i in range(0,len(DATA_Fig6a))],'kx')
    axarr[0].plot([DATA_Fig6a[i][0] for i in range(0,len(DATA_Fig6a))],[phosphs_nonnorm[0][i]/phosphs_nonnorm[0][0]*100 for i in range(0,len(DATA_Fig6a))],'k-')
    axarr[1].plot([DATA_Fig6b[i][0] for i in range(0,len(DATA_Fig6b))],[DATA_Fig6b[i][1] for i in range(0,len(DATA_Fig6b))],'kx')
    axarr[1].plot([DATA_Fig6b[i][0] for i in range(0,len(DATA_Fig6b))],[phosphs_nonnorm[1][i]/phosphs_nonnorm[1][0]*100 for i in range(0,len(DATA_Fig6b))],'k-')

    for ikey in range(0,len(VARIABLES)):
      VARmeaning = ['NgPKC_coeff','NgCaM_PKC_cat_coeff','Ng_dephosph'][ikey]
    #  axarr[0].text(0.75*DATA_Fig6a[-1][0], (0.9-0.08*ikey)*max([DATA_Fig6a[i][1] for i in range(0,len(DATA_Fig6a))]), VARIABLES[ikey][0]+'='+'{0:.7f}'.format(parameters[VARIABLES[ikey][0]])+'. '+VARmeaning,fontsize=5)
    #  axarr[1].text(0.75*DATA_Fig6b[-1][0], (0.9-0.08*ikey)*max([DATA_Fig6b[i][1] for i in range(0,len(DATA_Fig6b))]), '{0:.7f}'.format(10**parameters[VARIABLES[ikey][0]]),fontsize=5)
    #axarr[0].text(0.75*DATA_Fig6a[-1][0], 0.2*max([DATA_Fig6a[i][1] for i in range(0,len(DATA_Fig6a))]), OBJECTIVES[0]+'='+'{0:.7f}'.format(mydict[OBJECTIVES[0]]),fontsize=5,color='#FF0000')
    #axarr[1].text(0.75*DATA_Fig6a[-1][0], 0.2*max([DATA_Fig6b[i][1] for i in range(0,len(DATA_Fig6b))]), OBJECTIVES[1]+'='+'{0:.7f}'.format(mydict[OBJECTIVES[1]]),fontsize=5,color='#FF0000')
    Ngps = [sum([lastConcs_all[j][i] for i in range(0,len(headers)) if 'Ngp' in headers[i]]) for j in range(0,len(lastConcs_all))]
    Ngs = [sum([lastConcs_all[j][i] for i in range(0,len(headers)) if 'Ng' in headers[i]]) for j in range(0,len(lastConcs_all))]
    NgCaMs = [sum([lastConcs_all[j][i] for i in range(0,len(headers)) if 'Ng' in headers[i] and 'CaM' in headers[i]]) for j in range(0,len(lastConcs_all))]
    NgCaMCas = [sum([lastConcs_all[j][i] for i in range(0,len(headers)) if 'Ng' in headers[i] and 'CaMCa' in headers[i]]) for j in range(0,len(lastConcs_all))]
    CaMs = [sum([lastConcs_all[j][i] for i in range(0,len(headers)) if 'CaM' in headers[i]]) for j in range(0,len(lastConcs_all))]
    CaMCas = [sum([lastConcs_all[j][i] for i in range(0,len(headers)) if 'CaMCa' in headers[i]]) for j in range(0,len(lastConcs_all))]

    axarr[0].set_title('Normal Ca$^{2+}$',fontsize=6)
    axarr[1].set_title('With EGTA2', fontsize=6)
    axarr[0].set_xlabel('[CaM] (mM)', fontsize=6)
    axarr[1].set_xlabel('[CaM] (mM)', fontsize=6)
    axarr[0].set_ylabel('%', fontsize=6, rotation=0)
    axarr[1].set_ylabel('%', fontsize=6, rotation=0)

    f.savefig(saveFig+'.eps')
  return mydict
  
# After each generation this function is called
def checkpopulation(population, columns, gen):
  picklelist = [population, columns]
  print('Generation '+str(gen)+' done, saving to fitfiles/EGTA2_alt_seed'+str(myseed)+'_N'+str(N_samples)+'_tmp'+str(gen)+'.sav')
  file=open('fitfiles/EGTA2_alt_seed'+str(myseed)+'_N'+str(N_samples)+'_tmp'+str(gen)+'.sav', 'wb')
  pickle.dump(picklelist,file)
  file.close() 

C_samples = 2*N_samples
N_generations = 25

random.seed(19288+myseed)

summeds = [[i] for i in range(0,len(OBJECTIVES))] + [[i,j] for i,j in array(meshgrid(*[list(range(0,len(OBJECTIVES))),list(range(0,len(OBJECTIVES)))])).T.reshape(len(OBJECTIVES)**2,2) if i<j] +\
          [[i,j,k] for i,j,k in array(meshgrid(*[list(range(0,len(OBJECTIVES))),list(range(0,len(OBJECTIVES))),list(range(0,len(OBJECTIVES)))])).T.reshape(len(OBJECTIVES)**3,3) if i<j<k]

params_all = []
for gen in range(N_generations,-1,-1):
  if exists('fitfiles/EGTA2_alt_seed'+str(myseed)+'_N'+str(N_samples)+'_tmp'+str(gen)+'.sav'):
    unpicklefile = open('fitfiles/EGTA2_alt_seed'+str(myseed)+'_N'+str(N_samples)+'_tmp'+str(gen)+'.sav','rb');unpickledlist = pickle.load(unpicklefile,encoding='bytes');unpicklefile.close()
    params_all = params_all+unpickledlist[0].tolist()
  else:
    print('fitfiles/EGTA2_alt_seed'+str(myseed)+'_N'+str(N_samples)+'_tmp'+str(gen)+'.sav does not exist')
    

if True:
    params_all = array(params_all)
    param_fdims = list(range(len(VARIABLES), len(VARIABLES)+len(OBJECTIVES)))
    if unpickledlist[1][b'f0'] != len(VARIABLES):
      print("Something's wrong!"); time.sleep(3)
    medians0 = [median([y[i] for y in params_all]) for i in param_fdims]
    medians = [y if y > 0 else 1.0 for y in medians0]
    fvals = [sum([1.0*y[param_fdims[i]] for i in range(0,len(OBJECTIVES))]) for y in params_all]
    params_f =  [[1.0*y[param_fdims[j]]/medians[j] for j in range(0,len(param_fdims))] for y in params_all]

    for ifval in range(0,len(summeds)):
      if len(sys.argv) > 3 and int(sys.argv[3]) != ifval:
        continue
      myord = [i[0] for i in sorted(enumerate([sum([params_f[j][ifval2] for ifval2 in summeds[ifval]]) for j in range(0,len(params_f))]), key=lambda x:x[1])]
      for j in range(0,N_check):
        myparams={}; [myparams.update({'k'+str(i): params_all[myord[j]][i]}) for i in range(0,len(VARIABLES))]
        func_to_optimize(myparams,0.1,'figcaliD-E_seed'+str(myseed)+'_N'+str(N_samples)+'_'+chr(ord('a')+ifval)+str(j),0,True)
