#cp fitsheu_check_varEGTA_feweralt.py fitsheu_check_EGTA2.py
from pylab import *
import os
from os.path import exists
import random
import scipy.io
import emoo
import pickle
import mytools

Nstim = 100
freq = 100.0
Caflux = 200.0
Ntrains = 1
NgCoeff = 0.1


filenames = ['nrnhighres_tstop24640000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,Ngx0.5,0.5,1.5,1.5,'+str(NgCoeff)+'_k1x1.0_onset24040000.0_n'+str(Nstim)+'_freq'+str(freq)+'_dur3.0_flux'+str(Caflux)+'_Lflux5.0_Gluflux10.0_AChflux10.0_Ntrains'+str(Ntrains)+'_trainT4000.0.mat',
             'nrnhighres_tstop24640000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,Ngx0.5,0.5,1.5,1.5,'+str(NgCoeff)+'_k421,423x0.0,0.0_onset24040000.0_n'+str(Nstim)+'_freq'+str(freq)+'_dur3.0_flux'+str(Caflux)+'_Lflux5.0_Gluflux10.0_AChflux10.0_Ntrains'+str(Ntrains)+'_trainT4000.0.mat',
             'nrnhighres_tstop24640000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,Ngx0.5,0.5,1.5,1.5,'+str(NgCoeff)+'_k421,423x0.0,0.0_onset24040000.0_n'+str(Nstim)+'_freq'+str(freq)+'_dur3.0_flux'+str(Caflux)+'_Lflux5.0_Gluflux0.0_AChflux0.0_Ntrains'+str(Ntrains)+'_trainT4000.0.mat']
Glufluxes = [10.0, 10.0, 0.0]
#ks = ['1 1.0', '1 1.0', '421,423 0.0,0.0']
ks = ['1 1.0', '421,423 0.0,0.0', '421,423 0.0,0.0']
styles = ['k-','k--','k--']
mydashes = [(1,0),(4,1),(2,1)]

f,axarr = subplots(3,1)
for iax in range(0,3):
  axarr[iax].set_position([0.1+0.23*iax,0.8,0.15,0.15])
  axarr[iax].set_ylim([0,100])
  f.text(0.03+0.23*iax,0.94,chr(ord('I')+iax),fontsize=16)
  for tick in axarr[iax].xaxis.get_major_ticks() + axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(5)
  for axis in ['top','bottom','left','right']:
    axarr[iax].spines[axis].set_linewidth(0.5)

for ifile in range(0,len(filenames)):
  Gluflux = Glufluxes[ifile]
  myk = ks[ifile]
  if not exists(filenames[ifile]):
    print(filenames[ifile]+' does not exist')
    print('python3 model_nrn_altered_noU_highres.py 24640000 1e-6 24040000 '+str(Nstim)+' '+str(freq)+' 3.0 '+str(Caflux)+' 5.0 '+str(Gluflux)+' '+str(Gluflux)+' '+str(Ntrains)+' 4000 None GluR1,GluR1_memb,GluR2,GluR2_memb,Ng 0.5,0.5,1.5,1.5,'+str(NgCoeff)+' '+myk)
    os.system('python3 model_nrn_altered_noU_highres.py 24640000 1e-6 24040000 '+str(Nstim)+' '+str(freq)+' 3.0 '+str(Caflux)+' 5.0 '+str(Gluflux)+' '+str(Gluflux)+' '+str(Ntrains)+' 4000 None GluR1,GluR1_memb,GluR2,GluR2_memb,Ng 0.5,0.5,1.5,1.5,'+str(NgCoeff)+' '+myk)
  else:
    print(filenames[ifile]+' exists')

  A = scipy.io.loadmat(filenames[ifile])
  iPKCt_tot = [i for i in range(0,len(A['headers'])) if 'PKCt' in A['headers'][i]]
  iPKCt = [i for i in range(0,len(A['headers'])) if 'PKCt' in A['headers'][i] and 'Ng' not in A['headers'][i] and 'Glu' not in A['headers'][i]]
  #iPKCt = [i for i in range(0,len(A['headers'])) if A['headers'][i] == 'PKCtCa' or A['headers'][i] == 'PKCtAACa' or A['headers'][i] == 'PKCtDAGCa' or A['headers'][i] == 'PKCtAADAGCa']
  PKCt = sum(array([A['DATA'][i] for i in iPKCt]),axis=0)
  PKCt_tot = sum(array([A['DATA'][i] for i in iPKCt_tot]),axis=0)
  axarr[0].semilogx((A['DATA'][0]-24040000)/1000,PKCt_tot*1e6,styles[ifile],dashes=mydashes[ifile],lw=1.0)
  iabovehalf = [PKCt_tot[i] >= 0.5*max(PKCt_tot) for i in range(0,len(PKCt_tot))]
  firstbelow = [i for i in range(1,len(PKCt_tot)) if not iabovehalf[i] and iabovehalf[i-1]][0]
  firstabove = [i for i in range(1,len(PKCt_tot)) if iabovehalf[i] and not iabovehalf[i-1]][0]
  print('Half-activation time = '+str((A['DATA'][0][firstbelow]-24040000)/1000 - (A['DATA'][0][firstabove]-24040000)/1000))

  axarr[0].set_xlim([0,60]) #plot
  axarr[0].set_xlim([0.5,400]) #semilogx

axarr[0].set_title('[active PKC]',fontsize=6)
axarr[0].set_ylabel('nM',fontsize=6,rotation=0)
axarr[0].set_xlabel('time (s)',fontsize=6)
axarr[0].set_ylim([0,250]) #for Caflux 200.0
#axarr[0].set_ylim([0,2100]) #for Caflux 400.0

myleg = mytools.mylegend(f,[0.03,0.68,0.15,0.05],styles,['with DAG act., with AA bind.', 'with DAG act., no AA bind.', 'no DAG act., no AA bind.'],nx=1,dx=1.2,yplus=0.5,yplustext=0.35,colors=['#000000','#000000','#000000'],dashes=mydashes,linewidths=[1.0,1.0,1.0],myfontsize=6)
for q in ['top','right','bottom','left']:
  myleg.spines[q].set_visible(False)

myseed = 1
random.seed(71928+myseed)

LFlux = 1.0
dur = 10000.0



PDE4coeffs = [0.0, 0.5, 0.75, 1.0, 1.25, 1.50]
axarr[1].plot([1,6.0],[66.30564916637107/91.49151590398488*100]*2,'k--',color='#404040',lw=0.5)
MATs = []
PKAcs = []
xlabels = []
for iPDE4 in range(0,len(PDE4coeffs)):
  fixedBlock = 'PDE4'
  fixedBlockCoeff = str(PDE4coeffs[iPDE4])
  fixedAltered = '1'
  fixedAlteredCoeff = '1.0'
  filename = 'caliI_PDE4x'+str(PDE4coeffs[iPDE4])
  if not exists(filename+'.mat'):
    print("python3 model_nrn_altered_noU_extfilename.py 25000000 1e-6 24040000 1 1 10000.0 0.0 "+str(LFlux)+" 0.0 0.0 1 1 None "+fixedBlock+" "+fixedBlockCoeff+" "+fixedAltered+" "+fixedAlteredCoeff+" "+filename)
    os.system("python3 model_nrn_altered_noU_extfilename.py 25000000 1e-6 24040000 1 1 10000.0 0.0 "+str(LFlux)+" 0.0 0.0 1 1 None "+fixedBlock+" "+fixedBlockCoeff+" "+fixedAltered+" "+fixedAlteredCoeff+" "+filename)

  print("Loading "+filename+'.mat')
  MAT = scipy.io.loadmat(filename+'.mat')
  headers = MAT['headers']
  iPKAc = [i for i in range(0,len(headers)) if 'PKAc' == headers[i] or 'PKAc ' in headers[i]] #exclude 'PKAcAMP*', but include 'PKAc' with any number of whitespaces after 'c'
  PKAc = sum(array([MAT['DATA'][i,:] for i in iPKAc]),axis=0)
  PKAcs.append(PKAc[:])
  MATs.append(MAT)
  if iPDE4 > 0:
    mycolor = '#000000' if PDE4coeffs[iPDE4] == 1.0 else '#808080'              
    axarr[1].bar(0.5+iPDE4,max(PKAcs[-1])/max(PKAcs[0])*100,color = mycolor)
    xlabels.append([0.5+iPDE4,str(1*int((0.5+PDE4coeffs[iPDE4]*2814)/1))+' nM'])

axarr[1].set_title('PKA activity in WT\nrelative to PDE4-blocked', fontsize=6)
axarr[1].set_ylabel('%', fontsize=6, rotation=0)
axarr[1].set_xticks([xlabels[i][0] for i in range(0,len(xlabels))])
axarr[1].set_xticklabels([xlabels[i][1] for i in range(0,len(xlabels))],fontsize=6,rotation=30,ha='right')
axarr[1].set_xlabel('[PDE4] (nM)', fontsize=6)

f.savefig("figcaliI-K.eps")


PPcoeffs = [0.8, 0.9, 1.0, 1.1, 1.2]
axarr[2].plot([0,5.0],[15]*2,'k--',color='#404040',lw=0.5)
S845s = []
xlabels = []
for iPP in range(0,len(PPcoeffs)):
  fixedBlock = 'Ca,PP1,PP2B'
  fixedBlockCoeff = '1.0,'+str(PPcoeffs[iPP])+','+str(PPcoeffs[iPP])
  fixedAltered = '1'
  fixedAlteredCoeff = '1.0'
  filename = 'caliI_PPcoeffx'+str(PPcoeffs[iPP])
  if not exists(filename+'.mat'):
    print("python3 model_nrn_altered_noU_extfilename.py 25000000 1e-6 24040000 1 1 10000.0 0.0 "+str(LFlux)+" 0.0 0.0 1 1 None "+fixedBlock+" "+fixedBlockCoeff+" "+fixedAltered+" "+fixedAlteredCoeff+" "+filename)
    os.system("python3 model_nrn_altered_noU_extfilename.py 25000000 1e-6 24040000 1 1 10000.0 0.0 "+str(LFlux)+" 0.0 0.0 1 1 None "+fixedBlock+" "+fixedBlockCoeff+" "+fixedAltered+" "+fixedAlteredCoeff+" "+filename)

  print("Loading "+filename+'.mat')
  MAT = scipy.io.loadmat(filename+'.mat')
  headers = MAT['headers']
  iS845 = [i for i in range(0,len(headers)) if 'S845' in headers[i]] #exclude 'PKAcAMP*', but include 'PKAc' with any number of whitespaces after 'c'
  iGluR1 = [i for i in range(0,len(headers)) if 'GluR1' in headers[i]] #exclude 'PKAcAMP*', but include 'PKAc' with any number of whitespaces after 'c'
  S845 = sum(array([MAT['DATA'][i,:] for i in iS845]),axis=0)
  GluR1 = sum(array([MAT['DATA'][i,:] for i in iGluR1]),axis=0)
  S845s.append(S845[:]/GluR1[0])
  MATs.append(MAT)

  mycolor = '#000000' if PPcoeffs[iPP] == 1.0 else '#808080'              
  axarr[2].bar(0.5+iPP,S845[0]/GluR1[0]*100,color = mycolor)

  xlabels.append([0.5+iPP,str(int(PPcoeffs[iPP]*470+0.99))+'; '+str(int(PPcoeffs[iPP]*1080+0.99))])
  #axarr[2].text(0.5+iPP,10,str(PPcoeffs[iPP]*1.43e-08).replace('00000000000002',''),ha='center',va='bottom',color='#FFFFFF',rotation=90,fontsize=5)

axarr[2].set_title('Basal S845 phos.', fontsize=6)
axarr[2].set_ylabel('%', fontsize=6, rotation=0)
axarr[2].set_xticks([])
axarr[2].set_xticks([xlabels[i][0] for i in range(0,len(xlabels))])
axarr[2].set_xticklabels([xlabels[i][1]+' nM' for i in range(0,len(xlabels))],fontsize=6,rotation=30,ha='right')
axarr[2].set_xlabel('[PP1]; [PP2B] (nM)', fontsize=6)


#This won't be done!!!
#membbindcoeffs = [0.1, 0.3, 1.0, 3.0, 10.0]
#axarr[2].plot([0,5.0],[30]*2,'k--',color='#404040',lw=0.5)
#GluR1_membs = []
#xlabels = []
#for imemb in range(0,len(membbindcoeffs)):
#  fixedBlock = 'Ca'
#  fixedBlockCoeff = '1.0'
#  fixedAltered = '330,331,332,333,334,335,336,337,338,339,340,341'
#  fixedAlteredCoeff = ','.join([str(membbindcoeffs[imemb]) for x in range(0,1+sum([1 for x in fixedAltered if x == ',']))])
#  filename = 'caliI_membbindcoeffx'+str(membbindcoeffs[imemb])
#  if not exists(filename+'.mat'):
#    print("python3 model_nrn_altered_noU_extfilename.py 25000000 1e-6 24040000 1 1 10000.0 0.0 "+str(LFlux)+" 0.0 0.0 1 1 None "+fixedBlock+" "+fixedBlockCoeff+" "+fixedAltered+" "+fixedAlteredCoeff+" "+filename)
#    os.system("python3 model_nrn_altered_noU_extfilename.py 25000000 1e-6 24040000 1 1 10000.0 0.0 "+str(LFlux)+" 0.0 0.0 1 1 None "+fixedBlock+" "+fixedBlockCoeff+" "+fixedAltered+" "+fixedAlteredCoeff+" "+filename)
#
#  print("Loading "+filename+'.mat')
#  MAT = scipy.io.loadmat(filename+'.mat')
#  headers = MAT['headers']
#  iGluR1_memb = [i for i in range(0,len(headers)) if 'GluR1_memb' in headers[i]] #exclude 'PKAcAMP*', but include 'PKAc' with any number of whitespaces after 'c'
#  iGluR1 = [i for i in range(0,len(headers)) if 'GluR1' in headers[i]] #exclude 'PKAcAMP*', but include 'PKAc' with any number of whitespaces after 'c'
#  GluR1_memb = sum(array([MAT['DATA'][i,:] for i in iGluR1_memb]),axis=0)
#  GluR1 = sum(array([MAT['DATA'][i,:] for i in iGluR1]),axis=0)
#  GluR1_membs.append(GluR1_memb[:]/GluR1[0])
#  MATs.append(MAT)
#
#  mycolor = '#000000' if membbindcoeffs[imemb] == 1.0 else '#808080'              
#  axarr[2].bar(0.5+imemb,GluR1_memb[-1]/GluR1[0]*100,color = mycolor)
#  xlabels.append([0.5+imemb,str(membbindcoeffs[imemb]*0.95e-06).replace('00000000000005','').replace('00000000000002','').replace('499999999999997','5').replace('e-05','$\\times10^{-5}$').replace('e-06','$\\times10^{-6}$').replace('e-07','$\\times10^{-7}$')])
#  #axarr[2].text(0.5+imemb,7,str(membbindcoeffs[imemb]*0.95e-06).replace('00000000000005','').replace('00000000000002',''),ha='center',va='bottom',color='#FFFFFF',rotation=90,fontsize=5)
#
#axarr[2].set_title('Basal GluR1 at memb.', fontsize=6)
#axarr[2].set_ylabel('%', fontsize=6, rotation=0)
#axarr[2].set_xticks([])
#axarr[2].set_xticks([xlabels[i][0] for i in range(0,len(xlabels))])
#axarr[2].set_xticklabels([xlabels[i][1] for i in range(0,len(xlabels))],fontsize=6,rotation=30,ha='right')
#axarr[2].set_xlabel('$k_f$ (nM$^{-1}$ms$^{-1}$)', fontsize=6)


f.savefig("figcaliI-K.eps")


  

  
