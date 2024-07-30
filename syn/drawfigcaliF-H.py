#copied from spineNg2p/fitseki_new_grid.py
#cp fitseki_altd_grid.py fitseki_altf_grid.py
#cp fitseki_altb.py fitseki_altc.py
#cp fitseki_alt.py fitseki_altb.py
#cp fitseki.py fitseki_alt.py
from pylab import *
import os
from os.path import exists
import random
import scipy.io
import emoo
import pickle
import mytools

#These give for some reason different [Ng] than in the paper, wonder why...
#daltonInug = 1.66054e-18
#conc_Ng = 40.0/(7600*daltonInug)*6.022e-23*1e6 # 40ug/ml -> 40.0/(7600*daltonInug) molecules/ml -> 40.0/(7600*daltonInug)*6.022e-23 mol/ml -> 40.0/(7600*daltonInug)*6.022e-23*1e6 mM
conc_Ng = 0.0042     #mM, from Seki et al. 1995 Fig 2 caption
conc_PKC = 0.0       #mM #0.0000104 #mM, from Seki et al. 1995 (0.8 ug/l, calculated with molecular mass 77kDa) #Nope, PKC only in the phosphorylation assay, not anymore in the phosphatase assay!
conc_Ca = 0.5        #mM (TODO: check)
conc_CaM = 0.000588  #mM, from Seki et al. 1995 (10 ug/l, calculated with molecular mass 17kDa)
conc_PP1 = 30e-6/(37.5e3/6.022e23)/6.022e23*1e6 #mM (800nM):  g/ml -> N/ml -> mol/ml -> mM
conc_PP2A = 7.5e-6/(36e3/6.022e23)/6.022e23*1e6 #mM (208nM):  g/ml -> N/ml -> mol/ml -> mM
conc_PP2B = 4.7e-6/(59e3/6.022e23)/6.022e23*1e6 #mM (79.7nM): g/ml -> N/ml -> mol/ml -> mM
#conc_PP1 = 40e-6/(37.5e3/6.022e23)/6.022e23*1e6 #mM (800nM):  g/ml -> N/ml -> mol/ml -> mM
#conc_PP2A = 10e-6/(36e3/6.022e23)/6.022e23*1e6 #mM (208nM):  g/ml -> N/ml -> mol/ml -> mM
#conc_PP2B = 7.5e-6/(59e3/6.022e23)/6.022e23*1e6 #mM (79.7nM): g/ml -> N/ml -> mol/ml -> mM

DATA_Fig2a = [[5.435588802475152, 0.12359607111826887],  # PP1 dephos, Fig. 2A. Obtained with https://apps.automeris.io/wpd/
              [9.997033206891437, 0.21850961029053362],
              [19.47212274046918, 0.41190213821017596],
              [30.174193138231367, 0.6339613045413125],
              [59.05020237767277, 0.9389184979550318],
              [88.37186632477909, 0.9830442476000761],
              [118.76539023925068, 0.9946995062408611]]

DATA_Fig2b = [[5.684583718472712, 0.7469088733661715],   # PP2A, Fig. 2B
              [11.12643768342646, 0.9120631783867544],
              [20.514341895890993, 0.9928235354691919],
              [31.628128210346176, 0.9771735959064487],
              [60.00308928777206, 0.9836512929650054],
              [89.57551714922485, 0.9849483034661319],
              [117.78012607236289, 0.9879861031086572]]

DATA_Fig2c = [[5.020524341623364, 0.8220601031378891],   # PP2B, Fig. 2C
              [10.222154490529302, 0.87798217454785],
              [19.82636652787132, 0.921751869312791],
              [29.626727363003095, 0.9253848006803758],
              [59.367014556305506, 0.9345433272538677],
              [89.1167404887497, 0.9297403855135018],
              [119.69383089908904, 0.9511274868865434]]

fitsheu_alt_parameters = [-6.1665577,-5.3605618] # from check_EGTA2_alt_seed1_N1000_b1_obj.pdf

idata_seki = 0 #0 for 2a, 1 for 2b, 2 for 2c
myseed = 1
N_samples = 500

if len(sys.argv) > 1:
  idata_seki = int(float(sys.argv[1]))
if len(sys.argv) > 2:
  myseed = int(float(sys.argv[2]))
if len(sys.argv) > 3:
  N_samples = int(float(sys.argv[3]))

def phosphs(ks,idata_seki,verbose=False,rankID=0):
  NgpCaM_bind_coeff = 1.0
  NgpCaM_unbind_coeff = 1.0

  NgPKC_coeff = 10**fitsheu_alt_parameters[0]            #Multiply the backward rate of Ng + PKCt <-> NgPKCt by this coeff (and change forward and catalysis rates to keep the original Vmax and Km)
  NgCaM_PKC_bind_coeff = 1.0                             #Multiply the forward rate of NgCaM* + PKCt <-> NgPKCt by this coeff
  NgCaM_PKC_unbind_coeff = 1.0                           #Multiply the backward rate of NgCaM* + PKCt <-> NgPKCt by this coeff
  NgCaM_PKC_cat_coeff = 10**fitsheu_alt_parameters[1]    #Multiply the catalysis rate of NgCaM* + PKCt <-> NgPKCt by this coeff (compared to rate of Ng + PKCt <-> NgPKCt; Vmax and Km do not hold for NgCaM*)

  Ng_dephosph = 0.0

  Kforw = ks[0]
  Kbackw = ks[1]
  Kcat = ks[2]

  Ngps = []
  tvecs = []

  if idata_seki == 0:
    iKPP = '85,86,87'
  elif idata_seki == 1:
    iKPP = '88,89,90'
  elif idata_seki == 2:
    iKPP = '91,92,93'
  else:
    print("Unknown idata_seki")



  randomString = ''.join(random.choice('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789') for _ in range(0,10))+'_'+str(rankID)
  filenames = []

  filename = 'fit'+randomString+'.mat'

  #if idata_seki == 0:
  #  filename_IC = 'nrnICs_short_tstop0.0001_Ca-CaM-Ng-PKCt-PKCp-PKC-PP1-PP2A-PP2Bx'+str(1e6*conc_Ca)+'-'+str(1e6*conc_CaM)+'-'+str(1e6*conc_Ng)+'-'+str(0.5e6*conc_PKC)+'-'+str(0.5e6*conc_PKC)+'-0.0-0.0-'
  #
  #filename_IC = 'nrnICs_short_tstop0.0001_Ca-CaM-Ng-PKCt-PKCp-PKCx'+str(1e6*conc_Ca*iCaAmount)+'-'+str(1e6*conc_CaM)+'-'+str(1e6*conc_Ng)+'-'+str(0.5e6*conc_PKC)+'-'+str(0.5e6*conc_PKC)+'-0.0.mat'

  #Load the ICs from filename_IC. Remember that these will be in mM.
  #ICMAT = scipy.io.loadmat(filename_IC)
  #mystrIC_species = ''
  #mystrIC_concs = ''
  #for ispec2 in range(0,len(ICMAT['headers'])):
  #  header = ICMAT['headers'][ispec2]
  #  if header.find(' ') > -1:
  #    header = header[0:header.find(' ')]
  #  mystrIC_species = mystrIC_species + header + ','
  #  mystrIC_concs = mystrIC_concs + str(ICMAT['DATA'][ispec2][-1]) + ',' #These are printed in mM, like the data in ICMAT
  #mystrIC_species = mystrIC_species[0:-1]
  #mystrIC_concs = mystrIC_concs[0:-1]
  mystrIC_species = 'Ng,PKC,PKCp,Ngp,PKCt,Ca,CaM,PP1,PP2A,PP2B'
  mystrIC_concs = '0.0,0.0,0.0,'+str(conc_Ng)+','+str(conc_PKC)+','+str(conc_Ca)+','+str(conc_CaM)+','+str(conc_PP1 if idata_seki == 0 else 0.0)+','+str(conc_PP2A if idata_seki == 1 else 0.0)+','+str(conc_PP2B if idata_seki == 2 else 0.0)

  catperbackw = 887.2520420070011 # (kcat/kbackw = 1.7745040840140023 (1/ms) /0.002 (1/ms)
  NgPKC_coeff_forw = (NgPKC_coeff + catperbackw)/(1.0 + catperbackw)

  if idata_seki == 0:
    catperbackw_seki = 0.00035/0.0014 # (kcat/kbackw, unitless)
  elif idata_seki == 1:
    catperbackw_seki = 0.00015/0.005 # (kcat/kbackw, unitless)
  elif idata_seki == 2:
    catperbackw_seki = 0.002/0.008 # (kcat/kbackw, unitless)
  

  if verbose and not exists(filename):
    #      20,24,28,32,21,25,29,33,36,51,37,52,39,41,43,45,54,56,58,60,40,42,44,46,55,57,59,61,47,48,49,50,62,63,64,65,66,67,68,69,70
    print('python3 model_nrn_CaM_Ng_PKC_PPs_only_new_altered_extfilename_absconcs.py 8000000 1e-7 '+mystrIC_species+' '+mystrIC_concs+' '+\
          '34,38,42,46,35,39,43,47,50,65,51,66,53,55,57,59,68,70,72,74,54,56,58,60,69,71,73,75,61,62,63,64,76,77,78,79,80,81,82,83,84,'+iKPP+' '+\
           ','.join([str(NgpCaM_bind_coeff) for i in [0,1,2,3]])+','+','.join([str(NgpCaM_unbind_coeff) for i in [0,1,2,3]])+','+\
           ','.join([str(NgPKC_coeff_forw) for i in [0,1]])+','+','.join([str(NgPKC_coeff) for i in [0,1]])+','+\
           ','.join([str(NgCaM_PKC_bind_coeff*NgPKC_coeff_forw) for i in range(0,8)])+','+','.join([str(NgCaM_PKC_unbind_coeff*NgPKC_coeff) for i in range(0,8)])+','+','.join([str(NgCaM_PKC_cat_coeff) for i in range(0,8)])+','+\
           ','.join([str(Ng_dephosph) for i in range(0,5)])+','+str(Kforw)+','+str(Kbackw)+','+str(Kcat)+' '+filename)
    os.system('python3 model_nrn_CaM_Ng_PKC_PPs_only_new_altered_extfilename_absconcs.py 8000000 1e-7 '+mystrIC_species+' '+mystrIC_concs+' '+\
           '34,38,42,46,35,39,43,47,50,65,51,66,53,55,57,59,68,70,72,74,54,56,58,60,69,71,73,75,61,62,63,64,76,77,78,79,80,81,82,83,84,'+iKPP+' '+\
           ','.join([str(NgpCaM_bind_coeff) for i in [0,1,2,3]])+','+','.join([str(NgpCaM_unbind_coeff) for i in [0,1,2,3]])+','+\
           ','.join([str(NgPKC_coeff_forw) for i in [0,1]])+','+','.join([str(NgPKC_coeff) for i in [0,1]])+','+\
           ','.join([str(NgCaM_PKC_bind_coeff*NgPKC_coeff_forw) for i in range(0,8)])+','+','.join([str(NgCaM_PKC_unbind_coeff*NgPKC_coeff) for i in range(0,8)])+','+','.join([str(NgCaM_PKC_cat_coeff) for i in range(0,8)])+','+\
           ','.join([str(Ng_dephosph) for i in range(0,5)])+','+str(Kforw)+','+str(Kbackw)+','+str(Kcat)+' '+filename)
  elif not exists(filename):
    os.system('python3 model_nrn_CaM_Ng_PKC_PPs_only_new_altered_extfilename_absconcs.py 8000000 1e-7 '+mystrIC_species+' '+mystrIC_concs+' '+\
          '34,38,42,46,35,39,43,47,50,65,51,66,53,55,57,59,68,70,72,74,54,56,58,60,69,71,73,75,61,62,63,64,76,77,78,79,80,81,82,83,84,'+iKPP+' '+\
           ','.join([str(NgpCaM_bind_coeff) for i in [0,1,2,3]])+','+','.join([str(NgpCaM_unbind_coeff) for i in [0,1,2,3]])+','+\
           ','.join([str(NgPKC_coeff_forw) for i in [0,1]])+','+','.join([str(NgPKC_coeff) for i in [0,1]])+','+\
           ','.join([str(NgCaM_PKC_bind_coeff*NgPKC_coeff_forw) for i in range(0,8)])+','+','.join([str(NgCaM_PKC_unbind_coeff*NgPKC_coeff) for i in range(0,8)])+','+','.join([str(NgCaM_PKC_cat_coeff) for i in range(0,8)])+','+\
           ','.join([str(Ng_dephosph) for i in range(0,5)])+','+str(Kforw)+','+str(Kbackw)+','+str(Kcat)+' '+filename+' > /dev/null')
    ### These ks numbers should be added 14 !!! See difference between model_nrn_CaM_Ng_PKC_PPs_only_altered_extfilename_absconcs.py and model_nrn_CaM_Ng_PKC_only_altered_extfilename_absconcs.py ###
    #'20-24-28-32-' #Ngp binding to CaM*;     multiply by NgpCaM_bind_coeff
    #'21-25-29-33-' #Ngp unbinding from CaM*; multiply by NgpCaM_unbind_coeff; (this is alreaday tentatively 10 times faster than Ng unbinding from CaM*)
    #'36-51-' #Ng bind to PKC;     multiply by NgPKC_coeff_forw (ensures the right Km)
    #'37-52-' #Ng unbind from PKC; multiply by NgPKC_coeff
    #'38-53-' #NgPKC -> Ngp;       multiply by 1.0 (these rates are determined by Vmax and [PKCt], don't change them here)
    #'39-41-43-45-54-56-58-60-' #NgCaM bind to PKC;     multiply by NgCaM_PKC_bind_coeff*NgPKC_coeff_forw (correct Km if NgCaM_PKC_bind_coeff==1.0)
    #'40-42-44-46-55-57-59-61-' #NgCaM unbind from PKC; multiply by NgCaM_PKC_unbind_coeff*NgPKC_coeff    (correct Km if NgCaM_PKC_unbind_coeff==1.0)
    #'47-48-49-50-62-63-64-65-' #NgCaMPKC catalysis;    multiply by NgCaM_PKC_cat_coeff                   (correct Km if NgCaM_PKC_cat_coeff==1.0)
    #'66'                       # Ng_dephosph
    ### These ks numbers are OK !!! ###
    #'85'  # PP1 + Ngp forw
    #'86'  # PP1 + Ngp backw
    #'87'  # NgpPP1 cat
    #'88'  # PP2A + Ngp forw
    #'89'  # PP2A + Ngp backw
    #'90'  # NgpPP2A cat
    #'91'  # PP2B + Ngp forw
    #'92'  # PP2B + Ngp backw
    #'93'  # NgpPP2B cat
  if not exists(filename):
    Ngp = nan
    tvec = nan
  else:
    print("Loading "+filename)
    MAT = scipy.io.loadmat(filename)
    Ngp = sum(array([MAT['DATA'][j,:] for j in range(0,len(MAT['headers'])) if 'Ngp' in MAT['headers'][j]]),axis=0)
    tvec = MAT['DATA'][0]
    filenames.append(filename)

  #for filename in filenames:
  #  os.system('rm '+filename)

  return [Ngp, tvec]

OBJECTIVES = ['f0']
VARIABLES = [['k0',-5.,2.]], #Note: these are log10 values, i.e., each variable is allowed to vary from 0.0000001-fold to 100-fold
MAXERR = 1e8 # Maximum error for missing data

def func_to_optimize(params,idata_seki,saveFig="",rankID=0,verbose=False):
  
  coeffs = [10**params['k0'],1,1]
  if idata_seki == 0:
    DATA = array(DATA_Fig2a)
  elif idata_seki == 1:
    DATA = array(DATA_Fig2b)
  elif idata_seki == 2:
    DATA = array(DATA_Fig2c)

  phosphs_nonnorm_both = phosphs(coeffs, idata_seki, verbose, rankID=rankID)
  phosphs_nonnorm = phosphs_nonnorm_both[0] #the amount of phosphorylated Ng
  tvec = phosphs_nonnorm_both[1]/60000      #convert from ms to min
  try:
    DATA_interp = mytools.interpolate(tvec,phosphs_nonnorm,DATA[:,0])
  except:
    return {'f0': MAXERR}
  mydict = {'f0': sum((1-array(DATA_interp)/conc_Ng-DATA[:,1])**2)}
  if len(saveFig) > 0:
    f,axarr = subplots(1,1)
    axarr.plot(DATA[:,0],DATA[:,1],'bx')
    axarr.plot(tvec,1-phosphs_nonnorm/conc_Ng,'g-')
    axarr.set_title('f0 = '+str(mydict['f0']),fontsize=7)
    axarr.set_xlabel(str(params), fontsize=5)
    #axarr[1].set_title('f1 = '+str(mydict['f1']),fontsize=7)
    #axarr[1].set_xlabel('pars = '+str([10**parameters[VARIABLES[i][0]] for i in range(0,len(VARIABLES))]), fontsize=5)
    axarr.set_position([0.1,0.54,0.8,0.4])
    for tick in axarr.xaxis.get_major_ticks() + axarr.yaxis.get_major_ticks():
        tick.label.set_fontsize(4)
    f.savefig(saveFig+'_obj.eps')
  mydict['TARGETDATA'] = DATA
  mydict['phosphs'] = phosphs_nonnorm_both
  return mydict
  
# After each generation this function is called
def checkpopulation(population, columns, gen, idata_seki):
  picklelist = [population, columns]
  print('Generation '+str(gen)+' done, saving to fitfiles/sekinew'+str(idata_seki)+'_seed'+str(myseed)+'_N'+str(N_samples)+'_tmp'+str(gen)+'.sav')
  file=open('fitfiles/sekinew'+str(idata_seki)+'_seed'+str(myseed)+'_N'+str(N_samples)+'_tmp'+str(gen)+'.sav', 'w')
  pickle.dump(picklelist,file)
  file.close() 

C_samples = 2*N_samples
N_generations = 25

random.seed(23799 + 100*idata_seki + myseed)

data_all = []
f,axarr = subplots(3,1)
for iax in range(0,3):
  axarr[iax].set_position([0.1+0.23*iax,0.8,0.15,0.15])
  axarr[iax].set_ylim([0,100])
  f.text(0.03+0.23*iax,0.94,chr(ord('F')+iax),fontsize=16)
  for tick in axarr[iax].xaxis.get_major_ticks() + axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(5)
  for axis in ['top','bottom','left','right']:
    axarr[iax].spines[axis].set_linewidth(0.5)

mydicts = []
for idata in range(0,3):
  ks = [-3.0+0.1*i for i in range(0,50)]
  errs = []
  for ik in range(0,len(ks)):
    mydict = func_to_optimize({'k0': ks[ik]},idata,"test"+str(idata),0,True)
    errs.append(mydict['f0'])
  ibest = argmin(errs)
  ks2 = [-0.1+ks[ibest]+0.01*i for i in range(0,20)]
  errs2 = []
  for ik in range(0,len(ks2)):
    mydict = func_to_optimize({'k0': ks2[ik]},idata,"test"+str(idata),0,True)
    errs2.append(mydict['f0'])
  ibest2 = argmin(errs2)
  ks3 = [-0.01+ks2[ibest2]+0.001*i for i in range(0,20)]
  errs3 = []
  for ik in range(0,len(ks3)):
    mydict = func_to_optimize({'k0': ks3[ik]},idata,"test"+str(idata),0,True)
  errs3.append(mydict['f0'])
  ibest3 = argmin(errs3)
  k = ks3[ibest3]

  mydict = func_to_optimize({'k0': k},idata,"best_new"+str(idata),0,True)
  mydicts.append(mydict)
  data_all.append([ks,errs,ks2,errs2,ks3,errs3])

  DATA = mydict['TARGETDATA']
  phosphs_nonnorm_both = mydict['phosphs']
  phosphs_nonnorm = phosphs_nonnorm_both[0] #the amount of phosphorylated Ng
  tvec = phosphs_nonnorm_both[1]/60000      #convert from ms to min
  axarr[idata].plot(DATA[:,0],100*DATA[:,1],'kx')
  axarr[idata].plot(tvec,100*(1-phosphs_nonnorm/conc_Ng),'k-')
  #axarr.set_title('f0 = '+str(mydict['f0']),fontsize=7)
  #axarr.set_xlabel(str(params), fontsize=5)
  #axarr[1].set_title('f1 = '+str(mydict['f1']),fontsize=7)
  #axarr[1].set_xlabel('pars = '+str([10**parameters[VARIABLES[i][0]] for i in range(0,len(VARIABLES))]), fontsize=5)
  #axarr.set_position([0.1,0.54,0.8,0.4])
  axarr[idata].set_xlabel('time (min)', fontsize=6)
  axarr[idata].set_ylabel('%', fontsize=6, rotation=0)

axarr[0].set_title('PP1',fontsize=6)
axarr[1].set_title('PP2A',fontsize=6)
axarr[2].set_title('PP2B',fontsize=6)
f.savefig("figcaliF-H.eps")
