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

A = scipy.io.loadmat('nrn_tstop25840000_GluR1,GluR1_memb,GluR2,GluR2_memb,Cax0.5,0.5,1.5,1.5,1.0_n120_freq1.0_dur3.0_Lflux0.05_Gluflux0.0_AChflux0.05_Ntrains1_trainT100000.0_pair-100.0_icell1_pulseamp10.0_Nsyn1_Econ0.00015_wNMDA1.0_gNap0.03_Npulses4_apic250-300.mat') #Load this to find the time index
timevec = A['DATA'][0]
indT = argmin(abs(24040000+givenT*60000 - timevec))
print(str(givenT)+" min, indT = "+str(indT))

toBeRemovedIfNecessary = ['_CaCoeff1.0','_tol1e-06','_tstop15000000','3560000_600000','_Ninputs1','_tstop25000000','_onset24040000.0','_pulseamp5.0']

def boxoff(ax,whichxoff='top'):
    ax.spines[whichxoff].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def drawdiscontinuity(ax,y,yoffset,x=0,xoffset=0.1,lw=2.0,lw2=1.0,mycol='#000000'):
  thisline = ax.plot([x-xoffset,x+xoffset],[y-yoffset,y],'k-',linewidth=lw2,color=mycol)
  thisline[0].set_clip_on(False)
  thisline = ax.plot([x-xoffset,x+xoffset],[y,y+yoffset],'k-',linewidth=lw2,color=mycol)
  thisline[0].set_clip_on(False)
  thisline = ax.plot([x-1.2*xoffset,x+1.2*xoffset],[y-0.6*yoffset,y+0.6*yoffset],'k-',color='#FFFFFF',zorder=100,linewidth=lw)
  thisline[0].set_clip_on(False)

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

f,axs = subplots(5,3)
axarr = axs.reshape(prod(axs.shape),).tolist()
for iax in range(0,len(axarr)):
  boxoff(axarr[iax])
  if iax in [1,2,3,4,5,9,10]: #not in [0,6,7,8,11,12,13,14]: #here list of f-I plot indices
    axarr[iax].set_xlim([-90,170])
    axarr[iax].set_ylim([0.5,2.62])
    axarr[iax].plot([-250,250],[1.0,1.0],'k--',dashes=[4,2],lw=0.3)
  elif iax in [6,7,8]: #not in [11,12]:
    axarr[iax].set_xlim([0,0.25])
    axarr[iax].set_xticks([0,0.1,0.2])
    axarr[iax].set_ylim([0.0,5.0])
  elif iax in [11,12]:
    axarr[iax].set_ylim([0.0,3.0])
    axarr[iax].set_xlim([0,30])
  for tick in axarr[iax].xaxis.get_major_ticks()+axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(5)
  for axis in ['top','bottom','left','right']:
    axarr[iax].spines[axis].set_linewidth(0.3)
  axarr[iax].tick_params(axis='x',direction='out',width=0.4,length=2)
  axarr[iax].tick_params(axis='y',direction='out',width=0.4,length=2)

axarr[0].set_position([0.07,0.84,0.14,0.126])
axarr[0].set_xticks([])
axarr[0].set_yticks([])
boxoff(axarr[0])

for iax in range(0,2):
  axarr[1+iax].set_position([0.27+0.38*iax,0.81,0.32,0.126])
for iax in range(0,3):
  axarr[3+iax].set_position([0.07+0.32*iax,0.62,0.26,0.126])
  axarr[6+iax].set_position([0.07+0.32*iax,0.43,0.26,0.126])
for iax in range(0,2):
  axarr[9+iax].set_position([0.07+0.49*iax,0.24,0.42,0.126])
for iax in range(0,2):
  #axarr[11+iax].set_position([0.07+0.23*iax,0.81-0.19*4,0.17,0.14])
  axarr[11+iax].set_position([0.08+0.21*iax,0.05,0.19,0.126])
  axarr[13+iax].set_position([0.58+0.21*iax,0.05,0.19,0.126])
#axarr[21].set_position([0.75,0.81-0.25*3,0.2,0.14])

axinsets = []
for iax in range(0,15): #len(axarr)):
  axarr[iax].tick_params(axis='both', which='major', pad=2.0)
  axinsets.append(f.add_axes([axarr[iax].get_position().x0+0.045,axarr[iax].get_position().y0+0.065,0.03,0.056]))
  for tick in axinsets[iax].xaxis.get_major_ticks()+axinsets[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(5)
  for axis in ['top','bottom','left','right']:
    axinsets[iax].spines[axis].set_linewidth(0.3)
  if iax not in [6,7,8]: #insets for stdp curves
    boxoff(axinsets[iax]) #,'bottom')
    axinsets[iax].set_ylim([0,160])
    axinsets[iax].set_yticks([80,100,120]); axinsets[iax].set_yticklabels(['+80%','+100%','+120%']);
  else: #insets for f-I curves
    boxoff(axinsets[iax],'top')
    axinsets[iax].set_yticks([80,100,120]); axinsets[iax].set_yticklabels(['80%','100%','120%']);
  axinsets[iax].set_xticks([])
  axinsets[iax].tick_params(axis='y',direction='out',length=1,width=0.3)
  axinsets[iax].tick_params(axis='both', which='major', pad=0.1)
for iax in list(range(0,3))+list(range(11,15)):
  axinsets[iax].set_visible(False)
#for iax in [1,2,3,6,7,8,9,10,11,13,14]:
#  axarr[iax].set_yticklabels([])

########################## Draw the morphology of the L23PC #################################
#Pre-saved file with the data for plotting the morphologies
unpicklefile = open('morph_hilight250-300.sav', 'rb')
unpickledlist = pickle.load(unpicklefile,encoding='bytes')
unpicklefile.close()
segdata = unpickledlist[:]

for ipos in range(0,len(segdata)):
  coord1 = segdata[ipos][0]
  coord2 = segdata[ipos][1]
  mystyle = segdata[ipos][2]
  diam = segdata[ipos][3]
  mycol = segdata[ipos][4]
  if diam < 4:
    axarr[0].plot(coord1, coord2,mystyle,color=mycol, linewidth=diam*(0.5+0.5*(mycol != '#000000')))
  else:
    axarr[0].plot(coord1, coord2,mystyle,color=mycol, linewidth=1.0)
axarr[0].axis('equal')

polygon = Polygon(array([[0,295,270],[0,125,165]]).T, True)
p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
p.set_facecolor('#AAAAAA')
p.set_edgecolor('None')
axarr[0].add_collection(p)
polygon = Polygon(array([[-10,-310,-310],[250,227,273]]).T, True)
p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
p.set_facecolor('#AAAAAA')
p.set_edgecolor('None')
axarr[0].add_collection(p)
axarr[0].plot([-250,-250],[-100,0],'k-')
axarr[0].spines['left'].set_visible(False)
axarr[0].spines['bottom'].set_visible(False)
f.savefig('figstdp.eps')

conds,ISIs = get_stdp_curve('nrn_tstop25840000_GluR1,GluR1_memb,GluR2,GluR2_memb,Cax0.5,0.5,1.5,1.5,1.0_onset24040000.0_n120_freq1.0_dur3.0_Lflux0.0_Gluflux0.0_AChflux0.0_Ntrains1_trainT100000.0_pair0.0_icell1_pulseamp'+str(pulseamp)+'_Nsyn1_Econ'+str(Econ)+'_wNMDA1.0_gNap0.03_Npulses4_apic250-300.mat')
axarr[1].plot([x+30 for x in ISIs],conds,'b.-',mew=0.8,ms=0.8,lw=0.3,label='No $\\beta$-adr.')

conds,ISIs = get_stdp_curve('nrn_tstop25840000_GluR1,GluR1_memb,GluR2,GluR2_memb,Cax0.5,0.5,1.5,1.5,1.0_onset24040000.0_n120_freq1.0_dur3.0_Lflux0.05_Gluflux0.0_AChflux0.0_Ntrains1_trainT100000.0_pair0.0_icell1_pulseamp'+str(pulseamp)+'_Nsyn1_Econ'+str(Econ)+'_wNMDA1.0_gNap0.03_Npulses4_apic250-300.mat')
axarr[1].plot([x+30 for x in ISIs],conds,'k.-',mew=0.8,ms=0.8,lw=0.3,color='#964B00',label='With $\\beta$-adr.')

conds,ISIs = get_stdp_curve('nrn_tstop25840000_GluR1,GluR1_memb,GluR2,GluR2_memb,Cax0.5,0.5,1.5,1.5,1.0_onset24040000.0_n120_freq1.0_dur3.0_Lflux0.0_Gluflux0.0_AChflux0.05_Ntrains1_trainT100000.0_pair0.0_icell1_pulseamp'+str(pulseamp)+'_Nsyn1_Econ'+str(Econ)+'_wNMDA1.0_gNap0.03_Npulses4_apic250-300.mat')
axarr[2].plot([x+30 for x in ISIs],conds,'k.-',mew=0.8,ms=0.8,lw=0.3,color='#A020F0',label='No $\\beta$-adr.')

conds,ISIs = get_stdp_curve('nrn_tstop25840000_GluR1,GluR1_memb,GluR2,GluR2_memb,Cax0.5,0.5,1.5,1.5,1.0_onset24040000.0_n120_freq1.0_dur3.0_Lflux0.05_Gluflux0.0_AChflux0.05_Ntrains1_trainT100000.0_pair0.0_icell1_pulseamp'+str(pulseamp)+'_Nsyn1_Econ'+str(Econ)+'_wNMDA1.0_gNap0.03_Npulses4_apic250-300.mat')

axarr[2].plot([x+30 for x in ISIs],conds,'k.-',mew=0.8,ms=0.8,lw=0.3,label='With $\\beta$-adr.')
axarr[1].legend(fontsize=6,frameon=False)
axarr[2].legend(fontsize=6,frameon=False)
ISIs_control = ISIs[:]
conds_control = conds[:]
ts = [ISIs[i] for i in range(0,len(conds))]
tosumc = [conds_control[i] for i in range(0,len(conds_control)) if [x+30 for x in ISIs][i] >= -90 and [x+30 for x in ISIs][i] <= 170]
AUCc = sum([0.5*(tosumc[i]+tosumc[i+1])*(ts[i+1]-ts[i]) for i in range(0,len(tosumc)-1)])
spfreqs_control,amps_control = get_fi_curve('../l23pc/somaticDC_icell0_imut0_stimAmp0.0.mat')

cols = ['#666666', '#AA0000', '#FF00FF','#888800']
dimcols = ['#AAAAAA','#FF8888','#FF88FF','#AAAA44']
#  cols = ['#000000','#770077','#009900']
#  dimcols = ['#AAAAAA','#FF88FF','#99FF99']

blockedSet = [10,11]


conds,ISIs = get_stdp_curve('nrn_paired_imutCMcomb-1_comb10_'+str(Econ)+'_1.0_gNap0.03_pulseamp10.0_0.0.mat') #_imutCMcomb-1: no ion channel variants, _comb10: synaptic variants DLPFC, sig0, imp1, both
axarr[9].plot([x+30 for x in ISIs],conds,'k.-',mew=0.8,ms=0.8,lw=0.3,color=cols[2],label = 'PFC variant')
tosum = [conds[i] for i in range(0,len(conds)) if [x+30 for x in ISIs][i] >= -90 and [x+30 for x in ISIs][i] <= 170]
#axinsets[9].bar(1,(max(conds)-1)*100,color=cols[2]); print('PFC var baseline LTP '+str((max(conds)-1)*100))
axinsets[9].bar(1,sum([0.5*(tosum[i]+tosum[i+1])*(ts[i+1]-ts[i]) for i in range(0,len(tosum)-1)])/AUCc*100,color=cols[2])
if len(conds) > 0:
  print('PFC var baseline LTP '+str((max(conds)-1)*100))
  print("imutCMcomb-1, comb10 (no VGCC variants, syn PFC variants), max = "+str(max(conds))+", min = "+str(min(conds))+", conds[20] = "+str(conds[20])+", min(conds) = "+str(min(conds))+", imax = "+str(argmax(conds))+", imin = "+str(argmin(conds)))
conds,ISIs = get_stdp_curve('nrn_paired_imutCMcomb-1_comb11_'+str(Econ)+'_1.0_gNap0.03_pulseamp10.0_0.0.mat') #_imutCMcomb-1: no ion channel variants, _comb11: synaptic variants ACC, sig0, imp1, both
axarr[9].plot([x+30 for x in ISIs],conds,'k.-',mew=0.8,ms=0.8,lw=0.3,color=cols[3],label = 'ACC variant')
axarr[9].plot([x+30 for x in ISIs_control],conds_control,'k.-',mew=0.8,ms=0.8,lw=0.3,label = 'control')
tosum = [conds[i] for i in range(0,len(conds)) if [x+30 for x in ISIs][i] >= -90 and [x+30 for x in ISIs][i] <= 170]
if len(conds) > 0:
  print("imutCMcombr-1, comb11 (no VGCC variants, syn ACC variants), max = "+str(max(conds))+", min = "+str(min(conds))+", conds[20] = "+str(conds[20])+", min(conds) = "+str(min(conds))+", imax = "+str(argmax(conds))+", imin = "+str(argmin(conds)))
axinsets[9].bar(2,sum([0.5*(tosum[i]+tosum[i+1])*(ts[i+1]-ts[i]) for i in range(0,len(tosum)-1)])/AUCc*100,color=cols[3])
if len(conds) > 0:
  print('AnCg var baseline LTP '+str((max(conds)-1)*100))
axinsets[9].bar(0,100,color='#000000')
if len(conds_control) > 0:
  print('Control baseline LTP '+str((max(conds_control)-1)*100))

#Mutations (ion channel genes):
mutArray = ['gCa_HVAbar', 'gCa_LVAstbar', 'gIhbar', 'gK_Pstbar', 'gSK_E2bar',  'g', 'gNaTs2_tbar','gImbar','gK_Tstbar','gSKv3_1bar','decay', 'gamma', 'gNap_Et2bar']
mutSuffs = ['Ca_HVA', 'Ca_LVAst', 'Ih', 'K_Pst', 'SK_E2', 'pas', 'NaTs2_t', 'Im', 'K_Tst', 'SKv3_1', 'CaDynamics_E2', 'CaDynamics_E2', 'Nap_Et2']
titles = ['HVA Ca$^{2+}$', 'LVA Ca$^{2+}$', 'HCN', 'Slow K$^+$', 'SK', 'Leak', 'Fast Na$^+$', 'M-type K$^+$', 'Fast K$^+$', 'Kv3.1-type K$^+$', 'Ca$^{2+}$ decay', 'Ca$^{2+}$ fraction', 'Slow Na$^+$']
mutCoeffs = [0.75, 0.8, 0.9, 1.1, 1.2, 1.25]
#imutArray = [0, 1, 2, 3, 7, 6, 5]
imutArray = [0, 2]

data_for_corr = []
data_for_corr_max = []
for iimut in range(0,len(imutArray)):
  imut = imutArray[iimut]
  for imutcoeff in range(0,2):
    imutc = 6*imut+(2+3*imutcoeff)
    conds,ISIs = get_stdp_curve('nrn_tstop25840000_GluR1,GluR1_memb,GluR2,GluR2_memb,Cax0.5,0.5,1.5,1.5,1.0_onset24040000.0_n120_freq1.0_dur3.0_Lflux0.05_Gluflux0.0_AChflux0.05_Ntrains1_trainT100000.0_pair0.0_icell1_imutc'+str(imutc)+'_pulseamp'+str(pulseamp)+'_Nsyn1_Econ'+str(Econ)+'_wNMDA1.0_gNap0.03_Npulses4_apic250-300.mat')
    tosum = [conds[i] for i in range(0,len(conds)) if [x+30 for x in ISIs][i] >= -90 and [x+30 for x in ISIs][i] <= 170]
    spfreqs,amps = get_fi_curve('../l23pc/somaticDC_icell0_imutc'+str(imutc)+'_stimAmp0.0.mat')
    if len(conds) > 0:
      axarr[3+iimut].plot([x+30 for x in ISIs],conds,'k.-',mew=0.8,ms=0.8,color=cols[imutcoeff],lw=0.3, label = '-20%' if imutcoeff==1 else '+20%')
      axinsets[3+iimut].bar(imutcoeff*2,sum([0.5*(tosum[i]+tosum[i+1])*(ts[i+1]-ts[i]) for i in range(0,len(tosum)-1)])/AUCc*100,color=cols[imutcoeff])
    axarr[6+iimut].plot(amps,spfreqs,'k.-',mew=0.8,ms=0.8,color=cols[imutcoeff],lw=0.3, label = '-20%' if imutcoeff==1 else '+20%')
    if len(spfreqs) < 2:
      print('Something wrong with ../l23pc/somaticDC_icell0_imutc'+str(imutc)+'_stimAmp0.0.mat, skipping...')
      continue
    axinsets[6+iimut].bar(imutcoeff*2,sum(spfreqs)/sum(spfreqs_control)*100,color=cols[imutcoeff])
    tosum = [conds[i] for i in range(0,len(conds)) if [x+30 for x in ISIs][i] >= -90 and [x+30 for x in ISIs][i] <= 170]
    data_for_corr.append([sum([0.5*(tosum[i]+tosum[i+1])*(ts[i+1]-ts[i]) for i in range(0,len(tosum)-1)]),sum(spfreqs)*0.25]) #N.B. this only gives the correct correlation value when all variants are included, i.e. imutArray = [0, 1, 2, 3, 7, 6, 5]. Check drawfigstdp_withAUCs.py
    if len(conds) > 0 and len(spfreqs) > 0:
      data_for_corr_max.append([max(conds),max(spfreqs)])
    else:
      data_for_corr_max.append([nan,nan])

  axarr[3+iimut].plot([x+30 for x in ISIs_control],conds_control,'k.-',mew=0.8,ms=0.8,lw=0.3,label='control')
  axarr[6+iimut].plot(amps_control,spfreqs_control,'k.-',mew=0.8,ms=0.8,lw=0.3,label='control')
  axinsets[6+iimut].bar(1,100,color='#000000')
  axinsets[3+iimut].bar(1,100,color='#000000')

#Mutation combinations in ion channels, not in synaptic genes: STDP curve
conds,ISIs = get_stdp_curve('nrn_paired_imutCMcomb10_comb-1_'+str(Econ)+'_1.0_gNap0.03_pulseamp10.0_0.0.mat')
axarr[5].plot([x+30 for x in ISIs],conds,'k.-',mew=0.8,ms=0.8,color=cols[2],lw=0.3,label = 'PFC')
tosum = [conds[i] for i in range(0,len(conds)) if [x+30 for x in ISIs][i] >= -90 and [x+30 for x in ISIs][i] <= 170]
axinsets[5].bar(1,sum([0.5*(tosum[i]+tosum[i+1])*(ts[i+1]-ts[i]) for i in range(0,len(tosum)-1)])/AUCc*100,color=cols[2])
if len(conds) > 0:
  print("imutCMcomb10 (VGCC PFC variants), max = "+str(max(conds))+", min = "+str(min(conds))+", conds[20] = "+str(conds[20])+", min(conds) = "+str(min(conds))+", imax = "+str(argmax(conds))+", imin = "+str(argmin(conds)))
conds,ISIs = get_stdp_curve('nrn_paired_imutCMcomb11_comb-1_'+str(Econ)+'_1.0_gNap0.03_pulseamp10.0_0.0.mat')
axarr[5].plot([x+30 for x in ISIs],conds,'k.-',mew=0.8,ms=0.8,color=cols[3],lw=0.3,label = 'ACC')
tosum = [conds[i] for i in range(0,len(conds)) if [x+30 for x in ISIs][i] >= -90 and [x+30 for x in ISIs][i] <= 170]
axinsets[5].bar(2,sum([0.5*(tosum[i]+tosum[i+1])*(ts[i+1]-ts[i]) for i in range(0,len(tosum)-1)])/AUCc*100,color=cols[3])
axarr[5].plot([x+30 for x in ISIs_control],conds_control,'k.-',mew=0.8,ms=0.8,lw=0.3,label = 'control')
axinsets[5].bar(0,100,color='#000000')
if len(conds) > 0:
  print("imutCMcomb1 (VGCC ACC variants), max = "+str(max(conds))+", min = "+str(min(conds))+", conds[20] = "+str(conds[20])+", min(conds) = "+str(min(conds))+", imax = "+str(argmax(conds))+", imin = "+str(argmin(conds)))

#Mutation combinations in ion channels: f-I curves
spfreqs,amps = get_fi_curve('../l23pc/somaticDC_icell0_imutCMcomb10_stimAmp0.0.mat')
axarr[8].plot(amps,spfreqs,'k.-',mew=0.8,ms=0.8,color=cols[2],lw=0.3,label = 'PFC variant')
axinsets[8].bar(1,sum(spfreqs)/sum(spfreqs_control)*100,color=cols[2])
print("imutCMcomb10 (VGCC PFC variants), AUC = "+str(spfreqs)+" (len="+str(len(spfreqs))+"), perc = "+str(sum(spfreqs)/sum(spfreqs_control)*100))
spfreqs,amps = get_fi_curve('../l23pc/somaticDC_icell0_imutCMcomb11_stimAmp0.0.mat')
axarr[8].plot(amps,spfreqs,'k.-',mew=0.8,ms=0.8,color=cols[3],lw=0.3,label = 'PFC variant')
axinsets[8].bar(2,sum(spfreqs)/sum(spfreqs_control)*100,color=cols[3])
axarr[8].plot(amps_control,spfreqs_control,'k.-',mew=0.8,ms=0.8,lw=0.3,label = 'control')
axinsets[8].bar(0,100,color='#000000')
print("imutCMcomb11 (VGCC ACC variants), AUC = "+str(spfreqs)+" (len="+str(len(spfreqs))+"), perc = "+str(sum(spfreqs)/sum(spfreqs_control)*100))

#Mutation combinations in ion channels (DLPFC (combr0) and AnCg (combr1)), mutation combinations in synaptic genes (DLPFC (comb12) and AnCg (comb13))
conds,ISIs = get_stdp_curve('nrn_paired_imutCMcomb10_comb10_'+str(Econ)+'_1.0_gNap0.03_pulseamp10.0_0.0.mat') 
axarr[10].plot([x+30 for x in ISIs],conds,'k.-',mew=0.8,ms=0.8,color=cols[2],lw=0.3,label = 'PFC variant')
tosum = [conds[i] for i in range(0,len(conds)) if [x+30 for x in ISIs][i] >= -90 and [x+30 for x in ISIs][i] <= 170]
axinsets[10].bar(1,sum([0.5*(tosum[i]+tosum[i+1])*(ts[i+1]-ts[i]) for i in range(0,len(tosum)-1)])/AUCc*100,color=cols[2])
if len(conds) > 0:
  print("imutCMcomb10 (with VGCC variants), max = "+str(max(conds))+", min = "+str(min(conds))+", conds[20] = "+str(conds[20])+", min(conds) = "+str(min(conds))+", imax = "+str(argmax(conds))+", imin = "+str(argmin(conds)))
conds,ISIs = get_stdp_curve('nrn_paired_imutCMcomb11_comb11_'+str(Econ)+'_1.0_gNap0.03_pulseamp10.0_0.0.mat')
axarr[10].plot([x+30 for x in ISIs],conds,'k.-',mew=0.8,ms=0.8,color=cols[3],lw=0.3,label = 'ACC variant')
tosum = [conds[i] for i in range(0,len(conds)) if [x+30 for x in ISIs][i] >= -90 and [x+30 for x in ISIs][i] <= 170]
axinsets[10].bar(2,sum([0.5*(tosum[i]+tosum[i+1])*(ts[i+1]-ts[i]) for i in range(0,len(tosum)-1)])/AUCc*100,color=cols[3])
axarr[10].plot([x+30 for x in ISIs_control],conds_control,'k.-',mew=0.8,ms=0.8,lw=0.3,label = 'control')
axinsets[10].bar(0,100,color='#000000')
if len(conds) > 0:
  print("imutCMcomb11 (with VGCC variants), max = "+str(max(conds))+", min = "+str(min(conds))+", conds[20] = "+str(conds[20])+", min(conds) = "+str(min(conds))+", imax = "+str(argmax(conds))+", imin = "+str(argmin(conds)))

axarr[1].set_ylabel('Rel. cond. (A.U.)',fontsize=6.5)
axarr[3].set_ylabel('Rel. cond. (A.U.)',fontsize=6.5)
axarr[6].set_ylabel('firing rate (spikes/s)',fontsize=6.5)
#axarr[19].set_ylabel('Rel. cond. (A.U.)',fontsize=6.5)
#axarr[20].set_ylabel('firing rate (spikes/s)',fontsize=6.5)
#axarr[21].set_ylabel('Rel. cond. (A.U.)',fontsize=6.5)

for iax in [1,2,3,4,5,9,10]:
  axarr[iax].set_xlabel('ISI (ms)',fontsize=6.5,labelpad=1.0)
for iax in [6,7,8]:
  axarr[iax].set_xlabel('DC amplitude ($\mu$A)',fontsize=6.5,labelpad=1.0)
for iax in [11,12]:
  axarr[iax].set_xlabel('time (min)',fontsize=6.5,labelpad=1.0)
#axarr[19].set_xlabel('ISI (ms)',fontsize=6.5)
#axarr[20].set_xlabel('DC amplitude ($\mu$A)',fontsize=6.5)
#axarr[21].set_xlabel('ISI (ms)',fontsize=6.5)

axarr[1].set_title('No ACh',fontsize=6.5,pad=1.0)
axarr[2].set_title('With ACh',fontsize=6.5,pad=1.0)
for i in range(0,2):
  axarr[3+i].set_title(titles[imutArray[i]],fontsize=6.5,pad=1.0)
  axarr[6+i].set_title(titles[imutArray[i]],fontsize=6.5,pad=1.0)
axarr[5].set_title('Comb. (ion-channel encoding genes)',fontsize=6.5,pad=1.0)
axarr[8].set_title('Comb. (ion-channel encoding genes)',fontsize=6.5,pad=1.0)
axarr[9].set_title('Comb. (synaptic genes)',fontsize=6.5,pad=1.0)
axarr[10].set_title('Comb. (synaptic and ion-channel encoding genes)  ',fontsize=6.5,pad=1.0)

for iax in range(0,15):
  f.text(axarr[iax].get_position().x0-0.015,axarr[iax].get_position().y1+0.009-0.03*(iax==0),chr(ord('A')+iax),fontsize=8)



#print('Pearson corr between max LTP in STDP and max of fI curve: '+str(pearsonr(array(data_for_corr_max)[:,0],array(data_for_corr_max)[:,1])[0])+', pval = '+str(pearsonr(array(data_for_corr_max)[:,0],array(data_for_corr_max)[:,1])[1]))
#print('Pearson corr between max LTP in STDP and AUC of fI curve: '+str(pearsonr(array(data_for_corr_max)[:,0],array(data_for_corr)[:,1])[0])+', pval = '+str(pearsonr(array(data_for_corr_max)[:,0],array(data_for_corr)[:,1])[1]))
#print('Pearson corr between AUC LTP in STDP and max of fI curve: '+str(pearsonr(array(data_for_corr)[:,0],array(data_for_corr_max)[:,1])[0])+', pval = '+str(pearsonr(array(data_for_corr)[:,0],array(data_for_corr_max)[:,1])[1]))
#print('Pearson corr between AUC LTP in STDP and AUC of fI curve: '+str(pearsonr(array(data_for_corr)[:,0],array(data_for_corr)[:,1])[0])+', pval = '+str(pearsonr(array(data_for_corr)[:,0],array(data_for_corr)[:,1])[1]))
  
f.set_size_inches([6.4,6.4])
f.savefig('figstdp.eps')

### Sample-wise simulations

ISIs = [-10.0, -5.0, 0.0, 5.0, 10.0]
pulseamp = 10.0
gNap = 0.03
Econ = 0.00015

ISI_to_plot = 0.0

areas = ['PFC','ACC']
attrcommon = 'isSign0_isImputed1_Both'

for iax in range(0,2):
  axarr[11+iax].set_title(areas[iax]+", ISI="+str(ISI_to_plot+30)+" ms",fontsize=6.5,pad=1.0)
  axarr[13+iax].set_title(areas[iax],fontsize=6.5,pad=1.0)
axarr[11].set_ylabel('Relative synaptic\nstrength (fold change)',fontsize=6.5)
axarr[13].set_ylabel('Relative synaptic strength\nat 16 min (fold change)',fontsize=6.5)
axarr = axs.reshape(prod(axs.shape),).tolist()

def mybar(ax,x,y,facecolor=[],linewidth=0.3,w=0.4):
  qs = quantile(y, [0,0.25,0.5,0.75,1]) 
  polygon = Polygon(array([[x-w,x+w,x+w,x-w],[qs[1],qs[1],qs[3],qs[3]]]).T, True)
  p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
  if type(facecolor) is not list or len(facecolor) > 0:
    p.set_facecolor(facecolor)
  p.set_edgecolor('#000000')
  p.set_linewidth(0.3)
  ax.add_collection(p)

  a2 = ax.plot([x-w,x+w,x,x,x-w,x+w,x,x,x-w,x+w],[qs[0],qs[0],qs[0],qs[2],qs[2],qs[2],qs[2],qs[4],qs[4],qs[4]],'k-',lw=linewidth)
  return [p,a2]

for iarea in [0,1]:
  area = areas[iarea]
  attr = area+'_'+attrcommon
  if area == 'ACC':
    isubjs_HC = [1, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 73, 75, 76, 77, 78, 79, 80, 83, 84, 86, 90, 91, 92, 95, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 124, 126, 127, 129, 130, 131, 132, 133, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 152, 154, 159, 161, 163, 164, 167, 169, 170, 172, 174, 175, 176, 177, 182, 199, 208, 235, 240, 243, 249, 253, 256, 257, 260, 261, 262, 263, 265, 266, 267, 268, 271, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 329, 330, 331, 332, 343, 349, 350, 352, 353, 355, 360, 366, 369, 370, 371, 372, 374, 375, 376, 378, 381, 382, 383, 385, 388, 389, 392, 393, 395, 396, 397, 398, 401, 402, 409, 410, 414, 416, 418, 421, 422, 426, 427, 431, 432, 433, 434, 435, 437, 438, 442, 443, 444, 446, 447, 448, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480]
    isubjs_SCZ = [0, 2, 3, 12, 17, 28, 31, 58, 69, 70, 71, 72, 74, 81, 82, 85, 87, 88, 89, 93, 94, 96, 97, 98, 99, 100, 116, 117, 118, 119, 120, 121, 122, 123, 125, 134, 151, 155, 156, 157, 158, 160, 162, 165, 166, 168, 171, 173, 178, 179, 180, 181, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 200, 201, 202, 203, 204, 205, 206, 207, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 236, 237, 238, 239, 241, 242, 244, 245, 246, 247, 248, 250, 251, 252, 254, 255, 258, 259, 264, 269, 270, 272, 273, 274, 275, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 344, 345, 346, 347, 348, 351, 354, 356, 357, 358, 359, 361, 362, 363, 364, 365, 367, 368, 373, 377, 379, 380, 384, 386, 387, 390, 391, 394, 399, 400, 403, 404, 405, 406, 407, 411, 412, 413, 415, 417, 419, 420, 423, 424, 425, 428, 429, 430, 436, 439, 440, 441, 445, 449, 450, 451, 452]
  elif area == 'PFC':
    #From extract_commonmind_data_sklearn_allgenes_params_samples.py:
    isubjs_HC = [0, 1, 2, 3, 5, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 19, 20, 22, 23, 25, 27, 28, 31, 32, 33, 34, 35, 36, 37, 39, 41, 42, 43, 47, 51, 52, 53, 54, 63, 64, 66, 69, 73, 80, 81, 82, 83, 84, 87, 88, 89, 90, 92, 93, 94, 96, 97, 100, 105, 107, 111, 112, 118, 119, 120, 125, 127, 128, 135, 139, 140, 142, 143, 145, 146, 152, 154, 156, 165, 168, 176, 181, 182, 187, 192, 194, 201, 203, 204, 205, 207, 210, 211, 217, 220, 222, 223, 225, 226, 232, 235, 236, 240, 241, 242, 243, 244, 247, 248, 249, 251, 252, 253, 254, 257, 258, 262, 263, 265, 266, 267, 268, 269, 270, 271, 272, 283, 290, 291, 293, 294, 295, 296, 297, 298, 299, 300, 301, 305, 306, 308, 313, 314, 315, 317, 319, 320, 321, 323, 324, 325, 326, 329, 330, 332, 333, 335, 336, 339, 340, 341, 342, 343, 344, 345, 347, 351, 352, 354, 355, 356, 359, 361, 362, 368, 371, 373, 374, 375, 376, 377, 380, 381, 382, 383, 384, 386, 387, 388, 389, 390, 393, 396, 397, 398, 399, 400, 402, 403, 405, 406, 407, 408, 415, 416, 417, 420, 421, 422, 423, 424, 425, 426, 427, 428]
    isubjs_SCZ = [4, 6, 15, 17, 21, 24, 26, 29, 30, 38, 40, 44, 45, 46, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, 62, 65, 67, 68, 70, 71, 72, 74, 75, 76, 77, 78, 79, 85, 86, 91, 95, 98, 99, 101, 102, 103, 104, 106, 108, 109, 110, 113, 114, 115, 116, 117, 121, 122, 123, 124, 126, 129, 130, 131, 132, 133, 134, 136, 137, 138, 141, 144, 147, 148, 149, 150, 151, 153, 155, 157, 158, 159, 160, 161, 162, 163, 164, 166, 167, 169, 170, 171, 172, 173, 174, 175, 177, 178, 179, 180, 183, 184, 185, 186, 188, 189, 190, 191, 193, 195, 196, 197, 198, 199, 200, 202, 206, 208, 209, 212, 213, 214, 215, 216, 218, 219, 221, 224, 227, 228, 229, 230, 231, 233, 234, 237, 238, 239, 245, 246, 250, 255, 256, 259, 260, 261, 264, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 285, 286, 287, 288, 289, 292, 302, 303, 304, 307, 309, 310, 311, 312, 316, 318, 322, 327, 328, 334, 337, 338, 346, 348, 349, 350, 353, 357, 358, 360, 363, 364, 365, 366, 367, 369, 370, 372, 378, 379, 385, 392, 394, 395, 401, 404, 409, 410, 411, 412, 413, 414, 418, 419]
  else:
    error('Unknown area')
  #isubjs_HC = [0,1,2,3,4]; isubjs_SCZ = [5,6,7,8,9]
  
  condposts_thisarea = []
  condposts_abs_thisarea = []
  condtots_HC = []
  condtots_SCZ = []
  for iISI in range(0,len(ISIs)):
    ISI = ISIs[iISI]
    condpres = []
    condposts = []
    condposts_abs = []
    for isamp in range(0,481):
      filename =  'nrn_paired_imutCMcombcomb_'+area+'01Bothsamp'+str(isamp)+'_'+str(Econ)+'_1.0_gNap'+str(gNap)+'_pulseamp'+str(pulseamp)+'_'+str(ISI)+'.mat'
      if exists(filename):
        print('Loading '+filename)
        A = scipy.io.loadmat(filename)
        conds, times = calcconds.calcconds_nrn(filename)
        minabs = min(abs(times-(24040000+givenT*60*1000)))
        its = [it for it in range(0,len(times)) if times[it]-(24040000+givenT*60*1000) == minabs][0]
        condpres.append(conds[0])
        condposts.append(conds[its]/conds[0])
        condposts_abs.append(conds[its])
        if ISI == ISI_to_plot and isamp in isubjs_HC:
          condtots_HC.append(conds[:])
        if ISI == ISI_to_plot and isamp in isubjs_SCZ:
          condtots_SCZ.append(conds[:])
      else:
        print(filename+' does not exist')
        
    condposts_thisarea.append(condposts[:])
    condposts_abs_thisarea.append(condposts_abs[:])

    if len(condposts) == 0:
        continue

    pval = scipy.stats.ranksums(array([condposts[i] for i in isubjs_HC]),array([condposts[i] for i in isubjs_SCZ]))[1]
    print(area+',ISI='+str(ISI)+',av. condposts='+str(mean(array([condposts[i] for i in isubjs_HC])))+'+-'+str(std(array([condposts[i] for i in isubjs_HC])))+' vs '+str(mean(array([condposts[i] for i in isubjs_SCZ])))+'+-'+str(std(array([condposts[i] for i in isubjs_SCZ])))+', pval='+str(pval))

    mybar(axarr[13+iarea],3*iISI,[condposts[i] for i in isubjs_HC],facecolor='#666666')
    mybar(axarr[13+iarea],3*iISI+1,[condposts[i] for i in isubjs_SCZ],facecolor=cols[2+iarea])
    if pval < 0.05:
      #y1 = 1.03*max(mean([condposts[i] for i in isubjs_HC])+std([condposts[i] for i in isubjs_HC]),mean([condposts[i] for i in isubjs_SCZ])+std([condposts[i] for i in isubjs_SCZ]))
      y1 = 1.03*max(max([condposts[i] for i in isubjs_HC]),max([condposts[i] for i in isubjs_SCZ]))
      axarr[13+iarea].plot([3*iISI,3*iISI,3*iISI+1,3*iISI+1],[y1,y1*1.03,y1*1.03,y1],'k-',lw=0.23)
      axarr[13+iarea].plot(3*iISI+0.5,y1*1.08,'k*',mew=0.05,ms=3.0)
    #axarr[13+iarea].text(3*iISI+0.5,2.9+0.4*iarea,str(ISI+30)+" ms",ha='center',va='bottom',rotation=90,fontsize=6)
  axarr[13+iarea].set_xticks([3*iISI+0.5 for iISI in range(0,len(ISIs))])
  axarr[13+iarea].set_xticklabels([str(int(x+30)) for x in ISIs])
  axarr[13+iarea].set_xlabel('ISI (ms)',fontsize=6.5,labelpad=1.0)

  for isamp in range(0,len(condtots_HC)):
    axarr[11+iarea].plot((times-24040000)/60/1000,condtots_HC[isamp]/condtots_HC[isamp][0],color=dimcols[0],lw=0.23)
  for isamp in range(0,len(condtots_SCZ)):
    axarr[11+iarea].plot((times-24040000)/60/1000,condtots_SCZ[isamp]/condtots_SCZ[isamp][0],color=dimcols[2+iarea],lw=0.23)
  if len(condtots_HC) > 0:
    axarr[11+iarea].plot((times-24040000)/60/1000,[mean([condtots_HC[isamp][i]/condtots_HC[isamp][0] for isamp in range(0,len(condtots_HC))]) for i in range(0,len(condtots_HC[0]))],color='#000000',label='Control')
  if len(condtots_SCZ) > 0:
    axarr[11+iarea].plot((times-24040000)/60/1000,[mean([condtots_SCZ[isamp][i]/condtots_SCZ[isamp][0] for isamp in range(0,len(condtots_SCZ))]) for i in range(0,len(condtots_SCZ[0]))],color=cols[2+iarea],label='SCZ')
  axarr[11+iarea].legend(fontsize=6,loc=4,frameon=False,labelspacing=0.0,bbox_to_anchor=(1.05,-0.05))
  axarr[11+iarea].set_xlim([0,30])

for iax in [13,14]:
  axarr[iax].set_xlim([-0.75,13.75])
  axarr[iax].set_ylim([0.5,3.5])

for iax in range(0,len(axinsets)):
  ax = axinsets[iax]
  ax.set_ylim([73,125])
  if iax in [5,8,9,10]:
    cols = ['#000000','#FF00FF','#888800']
  else:
    cols = ['#666666', '#000000','#AA0000']

  for i in range(0,3):
    drawdiscontinuity(ax,78,2,i,0.4,lw=1.0,lw2=0.5,mycol=cols[i])
#axarr[5].legend(fontsize=6)
#myleg = mytools.mylegend(f,[0.025,0.77,0.2,0.044],['b-','b-','b-','b-','b-','b-'],['control','+20%','-20%','','PFC','ACC'],nx=2,dx=1.2,yplus=0.5,yplustext=0.35,colors=['#000000',cols[0],cols[1],'#FFFFFF',cols[2],cols[3]],myfontsize=6)
#myleg = mytools.mylegend(f,[0.025,0.77,0.25,0.02],['b-','b-','b-'],['control','+20%','-20%'],nx=3,dx=1.2,yplus=0.5,yplustext=0.35,colors=['#000000',cols[0],cols[1]],myfontsize=6)
#for q in ['top','bottom','left','right']:
#  myleg.spines[q].set_linewidth(0.2)
#myleg = mytools.mylegend(f,[0.81,0.625,0.1,0.02],['b-','b-'],['PFC','ACC'],nx=2,dx=1.2,yplus=0.5,yplustext=0.35,colors=[cols[2],cols[3]],myfontsize=6)
#for q in ['top','bottom','left','right']:
#  myleg.spines[q].set_linewidth(0.2)

handles, labels = axarr[3].get_legend_handles_labels()
if len(labels) > 2:
  leg = axarr[3].legend([handles[idx] for idx in [0,2,1]],[labels[idx] for idx in [0,2,1]],fontsize=6,frameon=False,loc=1,labelspacing=0.0,bbox_to_anchor=(1.1,1.0))
handles, labels = axarr[5].get_legend_handles_labels()
if len(labels) > 2:
  leg = axarr[5].legend([handles[idx] for idx in [2,0,1]],[labels[idx] for idx in [2,0,1]],fontsize=6,frameon=False,loc=1,labelspacing=0.0)
f.savefig("figstdp.eps")



