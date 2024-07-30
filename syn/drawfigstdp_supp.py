# cp ../spineM9l/drawfig6varcontnm.py drawfigstdp.py
#drawfig6.py: Draws the STDP curves
#Tuomo Maki-Marttunen, 2019-2020
#CC BY 4.0
import matplotlib
matplotlib.use('Agg')
from pylab import *
import scipy.io
import sys
import itertools
from os.path import exists
import mytools
import pickle
from matplotlib.collections import PatchCollection
import calcconds
from scipy.stats.stats import pearsonr

TRAINISISALL=[-230., -220., -210., -200., -180., -160., -140., -120., -100., -80., -60., -50., -40., -30., -25., -20., -15., -10., -5., -2.5, 0., 2.5, 5., 10., 15., 20., 25., 30., 40., 50., 60., 80., 100., 120., 140., 160., 180., 200. ]
Econ = 0.00015
pulseamp = 10.0
gNap = 0.03
V_min_peak=-40
V_max_valley=0
givenT = 16

A = scipy.io.loadmat('nrn_tstop25840000_GluR1,GluR1_memb,GluR2,GluR2_memb,Cax0.5,0.5,1.5,1.5,1.0_onset24040000.0_n120_freq1.0_dur3.0_Lflux0.05_Gluflux0.0_AChflux0.05_Ntrains1_trainT100000.0_pair0.0_icell1_imutc0_pulseamp5.0_Nsyn1_Econ0.00014_wNMDA1.0_Npulses4_apic250-300.mat') #Load this to find the time index
timevec = A['DATA'][0]
indT = argmin(abs(24040000+givenT*60000 - timevec))
print(str(givenT)+" min, indT = "+str(indT))

toBeRemovedIfNecessary = ['_CaCoeff1.0','_tol1e-06','_tstop15000000','3560000_600000','_Ninputs1','_tstop25000000','_onset24040000.0','_pulseamp5.0']

def boxoff(ax,whichxoff='top'):
    ax.spines[whichxoff].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def get_stdp_curve(filename):
  if filename.find('paired') > -1:
    filenames = [filename.replace('0.0.mat',str(x)+'.mat') for x in TRAINISISALL]
  else:
    filenames = [filename.replace('pair0.0','pair'+str(x)) for x in TRAINISISALL]
  ISIs = []
  conds = []
  for ifile in range(0,len(filenames)):
    if exists(filenames[ifile]):
      print('Loading '+filenames[ifile])
      cond, times = calcconds.calcconds_nrn(filenames[ifile])
      ISIs.append(TRAINISISALL[ifile])
      conds.append(cond[indT]/cond[0])
    else:
      myfilename = filenames[ifile]
      for i in range(0,len(toBeRemovedIfNecessary)):
        if len(myfilename) > 254:
          myfilename = myfilename.replace(toBeRemovedIfNecessary[i],'')
          #print("myfilename not short enough previously, now "+myfilename)
      if exists(myfilename):
        print('Loading '+myfilename)
        cond, times = calcconds.calcconds_nrn(myfilename)
        ISIs.append(TRAINISISALL[ifile])
        conds.append(cond[indT]/cond[0])
      else:
        print(filenames[ifile]+' does not exist')
  if len(ISIs) > 20 and ISIs[20] != 0.0:
    print ("Something wrong here")
  if len(conds) == 0:
    print(filename+": max cond not determined, len=0")
  else:
    print(filename+": max cond = "+str(max(conds)))
  return [conds, ISIs]

def get_fi_curve(filename):
  spfreqs = []
  amps_saved = []
  amps=[0.0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25]
  filenames = [filename.replace('stimAmp0.0','stimAmp'+str(x)) for x in amps]
  for iamp in range(0,len(amps)):
    if exists(filenames[iamp]):
      A = scipy.io.loadmat(filenames[iamp])
      sps = mytools.spike_times(A['times'][0],A['vsoma'][0],V_min_peak=V_min_peak,V_max_valley=V_max_valley)
      spfreqs.append(len(sps)/19.0) #[x/19.0 for x in Nsps]
      amps_saved.append(amps[iamp])
    else:
      print(filenames[iamp]+' does not exist')
  if len(amps_saved) > 0 and amps_saved[-1] != 0.25:
    print("Something wrong here, spfreqs")
  return [spfreqs, amps_saved]

f,axs = subplots(5,2)
axarr = sum([axs[i].tolist() for i in range(0,len(axs))]+[[]])
for iax in range(0,len(axarr)):
  boxoff(axarr[iax])
  if iax in [0,1,2,3,4]: 
    axarr[iax].set_position([0.08+0.19*iax,0.6,0.14,0.25])
    axarr[iax].set_xlim([-90,170])
    axarr[iax].set_ylim([0.5,2.62])
    axarr[iax].plot([-250,250],[1.0,1.0],'k--',dashes=[4,2],lw=0.3)
  elif iax in [5,6,7,8,9]:
    axarr[iax].set_position([0.08+0.19*(iax-5),0.27,0.14,0.25])
    axarr[iax].set_xlim([0,0.25])
    axarr[iax].set_xticks([0,0.1,0.2])
    axarr[iax].set_ylim([0.0,5.0])
  for tick in axarr[iax].xaxis.get_major_ticks()+axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(5)
  for axis in ['top','bottom','left','right']:
    axarr[iax].spines[axis].set_linewidth(0.3)
  axarr[iax].tick_params(axis='x',direction='out',width=0.4,length=2)
  axarr[iax].tick_params(axis='y',direction='out',width=0.4,length=2)

#for iax in range(0,2):
#  axarr[1+iax].set_position([0.27+0.38*iax,0.81,0.32,0.126])
#for iax in range(0,3):
#  axarr[3+iax].set_position([0.07+0.32*iax,0.62,0.26,0.126])
#  axarr[6+iax].set_position([0.07+0.32*iax,0.43,0.26,0.126])
#for iax in range(0,2):
#  axarr[9+iax].set_position([0.07+0.49*iax,0.24,0.42,0.126])
#for iax in range(0,2):
#  axarr[11+iax].set_position([0.08+0.21*iax,0.05,0.19,0.126])
#  axarr[13+iax].set_position([0.58+0.21*iax,0.01,0.19,0.166])

axinsets = []
for iax in range(0,len(axarr)):
  axarr[iax].tick_params(axis='both', which='major', pad=2.0)
  if iax < 5:
    axinsets.append(f.add_axes([axarr[iax].get_position().x0+0.13,axarr[iax].get_position().y0+0.17,0.02,0.056]))
  else:
    axinsets.append(f.add_axes([axarr[iax].get_position().x0+0.13,axarr[iax].get_position().y0+0.04,0.02,0.056]))
  for tick in axinsets[iax].xaxis.get_major_ticks()+axinsets[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(5)
  for axis in ['top','bottom','left','right']:
    axinsets[iax].spines[axis].set_linewidth(0.3)
  if iax not in [5,6,7,8,9]: #insets for stdp curves
    boxoff(axinsets[iax]) #,'bottom')
    axinsets[iax].set_ylim([0,160])
    axinsets[iax].set_yticks([50,100,150]); axinsets[iax].set_yticklabels(['+50%','+100%','+150%']);
  else: #insets for f-I curves
    boxoff(axinsets[iax],'top')
    axinsets[iax].set_ylim([0,125])
    axinsets[iax].set_yticks([0,50,100]); axinsets[iax].set_yticklabels(['0','50%','100%']);
  axinsets[iax].set_xticks([])
  axinsets[iax].tick_params(axis='y',direction='out',length=1,width=0.3)
  axinsets[iax].tick_params(axis='both', which='major', pad=0.1)
#for iax in list(range(0,3))+list(range(11,15)):
#  axinsets[iax].set_visible(False)


conds,ISIs = get_stdp_curve('nrn_tstop25840000_GluR1,GluR1_memb,GluR2,GluR2_memb,Cax0.5,0.5,1.5,1.5,1.0_onset24040000.0_n120_freq1.0_dur3.0_Lflux0.05_Gluflux0.0_AChflux0.05_Ntrains1_trainT100000.0_pair0.0_icell1_pulseamp'+str(pulseamp)+'_Nsyn1_Econ'+str(Econ)+'_wNMDA1.0_gNap0.03_Npulses4_apic250-300.mat')

ISIs_control = ISIs[:]
conds_control = conds[:]
spfreqs_control,amps_control = get_fi_curve('../l23pc4/somaticDC_icell0_imut0_stimAmp0.0.mat')
print("control, max = "+str(max(conds))+", min = "+str(min(conds))+", conds[20] = "+str(conds[20])+", min(conds) = "+str(min(conds))+", imax = "+str(argmax(conds))+", imin = "+str(argmin(conds)))

cols = ['#666666', '#AA0000', '#FF00FF','#888800']
dimcols = ['#AAAAAA','#FF8888','#FF88FF','#AAAA44']

blockedSet = [10,11]

#Mutations (ion channel genes):
mutArray = ['gCa_HVAbar', 'gCa_LVAstbar', 'gIhbar', 'gK_Pstbar', 'gSK_E2bar',  'g', 'gNaTs2_tbar','gImbar','gK_Tstbar','gSKv3_1bar','decay', 'gamma', 'gNap_Et2bar']
mutSuffs = ['Ca_HVA', 'Ca_LVAst', 'Ih', 'K_Pst', 'SK_E2', 'pas', 'NaTs2_t', 'Im', 'K_Tst', 'SKv3_1', 'CaDynamics_E2', 'CaDynamics_E2', 'Nap_Et2']
titles = ['HVA Ca$^{2+}$', 'LVA Ca$^{2+}$', 'HCN', 'Slow K$^+$', 'SK', 'Leak', 'Fast Na$^+$', 'M-type K$^+$', 'Fast K$^+$', 'Kv3.1-type K$^+$', 'Ca$^{2+}$ decay', 'Ca$^{2+}$ fraction', 'Slow Na$^+$']
mutCoeffs = [0.75, 0.8, 0.9, 1.1, 1.2, 1.25]
#imutArray = [0, 1, 2, 3, 7, 6, 5]
imutArray = [7, 5, 6, 3, 1]

data_for_corr = []
data_for_corr_max = []
for iimut in range(0,len(imutArray)):
  imut = imutArray[iimut]
  for imutcoeff in range(0,2):
    imutc = 6*imut+(2+3*imutcoeff)
    conds,ISIs = get_stdp_curve('nrn_tstop25840000_GluR1,GluR1_memb,GluR2,GluR2_memb,Cax0.5,0.5,1.5,1.5,1.0_onset24040000.0_n120_freq1.0_dur3.0_Lflux0.05_Gluflux0.0_AChflux0.05_Ntrains1_trainT100000.0_pair0.0_icell1_imutc'+str(imutc)+'_pulseamp'+str(pulseamp)+'_Nsyn1_Econ'+str(Econ)+'_wNMDA1.0_gNap0.03_Npulses4_apic250-300.mat')
    spfreqs,amps = get_fi_curve('../l23pc4/somaticDC_icell0_imutc'+str(imutc)+'_stimAmp0.0.mat')
    if len(conds) > 0:
      axarr[iimut].plot([x+30 for x in ISIs],conds,'k.-',mew=0.8,ms=0.8,color=cols[imutcoeff],lw=0.3, label = '-20%' if imutcoeff==1 else '+20%')
      axinsets[iimut].bar(imutcoeff*2,(max(conds)-1)*100,color=cols[imutcoeff])
    axarr[5+iimut].plot(amps,spfreqs,'k.-',mew=0.8,ms=0.8,color=cols[imutcoeff],lw=0.3, label = '-20%' if imutcoeff==1 else '+20%')
    if len(spfreqs) < 2:
      print('Something wrong with ../l23pc4/somaticDC_icell0_imutc'+str(imutc)+'_stimAmp0.0.mat, skipping...')
      continue
    axinsets[5+iimut].bar(imutcoeff*2,spfreqs[-1]/spfreqs_control[-1]*100,color=cols[imutcoeff])
    tosum = [conds[i] for i in range(0,len(conds)) if [x+30 for x in ISIs][i] >= -90 and [x+30 for x in ISIs][i] <= 170]
    ts = [ISIs[i] for i in range(0,len(conds)) if [x+30 for x in ISIs][i] >= -90 and [x+30 for x in ISIs][i] <= 170]    
    data_for_corr.append([sum([0.5*(tosum[i]+tosum[i+1])*(ts[i+1]-ts[i]) for i in range(0,len(tosum)-1)]),sum(spfreqs)*0.25])
    data_for_corr_max.append([max(conds),max(spfreqs)*0.25])

  axarr[iimut].plot([x+30 for x in ISIs_control],conds_control,'k.-',mew=0.8,ms=0.8,lw=0.3,label='control')
  axarr[5+iimut].plot(amps_control,spfreqs_control,'k.-',mew=0.8,ms=0.8,lw=0.3,label='control')
  axinsets[5+iimut].bar(1,100,color='#000000')
  axinsets[iimut].bar(1,(max(conds_control)-1)*100,color='#000000')


axarr[0].set_ylabel('Rel. cond. (A.U.)',fontsize=6.5)
axarr[5].set_ylabel('firing rate (spikes/s)',fontsize=6.5)

for iax in [0,1,2,3,4]:
  axarr[iax].set_xlabel('ISI (ms)',fontsize=6.5,labelpad=1.0)
for iax in [5,6,7,8,9]:
  axarr[iax].set_xlabel('DC amplitude ($\mu$A)',fontsize=6.5,labelpad=1.0)

for i in range(0,5):
  axarr[i].set_title(titles[imutArray[i]],fontsize=6.5,pad=1.0)
  axarr[5+i].set_title(titles[imutArray[i]],fontsize=6.5,pad=1.0)

for iax in range(0,10):
  f.text(axarr[iax].get_position().x0-0.015,axarr[iax].get_position().y1+0.009-0.03*(iax==0),chr(ord('A')+iax),fontsize=8)

f.savefig('figstdp_supp.eps')




