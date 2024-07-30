#cp drawfig_ltdltpcurves_muts_all_ramakercomb.py drawfig_ltdltpcurves_muts_all.py #7.11.2022
from pylab import *
import scipy.io
from os.path import exists
import os
import mytools
import calcconds
from matplotlib.collections import PatchCollection
import scipy.stats

#boxoff: remove extra axes (top and right) from an axis
def boxoff(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

#ltdltpcurve: retrieve data from completed simulations for a single variant
def ltdltpcurve(blocked='GluR1,GluR1_memb,GluR2,GluR2_memb,Cax0.5,0.5,1.5,1.5,1.0', T = 100, Caflux = 100.0, Lflux = 5.0, Gluflux = 10.0, altered = '_k1x1.0', freqs = [0.5, 1.0, 2.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 22.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 200.0, 300.0, 500.0]):

 Nstims = [int(T*x) for x in freqs]

 tposts = [30, 0] #min after stimulus onset

 condposts = []
 condposts_abs = []
 GluR1posts = []
 GluR1S831posts = []
 GluR2posts = []


 addition = ''
 for ifreq in range(0,len(freqs)):
  freq = freqs[ifreq]
  Nstim = Nstims[ifreq]
  filename = addition+'nrn_tstop27000000_tol1e-06_'+blocked+altered+'_onset24040000.0_n'+str(Nstim)+'_freq'+str(freq)+'_dur3.0_flux'+str(Caflux)+'_Lflux'+str(Lflux)+'_Gluflux'+str(Gluflux)+'_AChflux'+str(Gluflux)+'_Ntrains1_trainT1000.0.mat'
  if ifreq == 0 and not exists(filename) and not exists(filename+'.mat') and (exists('muts/'+filename) or exists('muts/'+filename+'.mat')): #do this only once to save i/o, the same should apply to either all freqs or none
    addition = 'muts/'
    filename = 'muts/'+filename

  if exists(filename):
    print('Loading '+filename)
    conds, times = calcconds.calcconds_nrn(filename)
    minabs = [min(abs(times-(24040000+x*60*1000))) for x in tposts]
    its = [[it for it in range(0,len(times)) if times[it]-(24040000+tposts[ipost]*60*1000) == minabs[ipost]][0] for ipost in range(0,len(tposts))]
    condposts.append([conds[its[ipost]]/conds[0] for ipost in range(0,len(tposts))])
    condposts_abs.append([conds[its[ipost]] for ipost in range(0,len(tposts))])
    A = scipy.io.loadmat(filename)
    iGluR1 = [i for i in range(0,len(A['headers'])) if 'GluR1_memb' in A['headers'][i]]
    iGluR1S831 = [i for i in range(0,len(A['headers'])) if 'GluR1_memb' in A['headers'][i] and 'S831' in A['headers'][i]]
    iGluR2 = [i for i in range(0,len(A['headers'])) if 'GluR2_memb' in A['headers'][i]]
    GluR1s = [sum([A['DATA'][i][it] for i in iGluR1]) for it in range(0,len(times))]
    GluR1S831s = [sum([A['DATA'][i][it] for i in iGluR1S831]) for it in range(0,len(times))]
    GluR2s = [sum([A['DATA'][i][it] for i in iGluR2]) for it in range(0,len(times))]
    GluR1posts.append([GluR1s[its[ipost]] for ipost in range(0,len(tposts))])
    GluR1S831posts.append([GluR1S831s[its[ipost]] for ipost in range(0,len(tposts))])
    GluR2posts.append([GluR2s[its[ipost]] for ipost in range(0,len(tposts))])
  else:
    print(filename+' not found')
    condposts.append(list(nan*ones([len(tposts),])))
    condposts_abs.append(list(nan*ones([len(tposts),])))
    GluR1posts.append(list(nan*ones([len(tposts),])))
    GluR1S831posts.append(list(nan*ones([len(tposts),])))
    GluR2posts.append(list(nan*ones([len(tposts),])))

 return [condposts, condposts_abs, GluR1posts, GluR2posts, GluR1S831posts]
  
#ltdltpcurve: retrieve data from completed simulations for a group of variants
def ltdltpcurves_many(blockeds=['GluR1,GluR1_memb,GluR2,GluR2_memb,Cax0.5,0.5,1.5,1.5,1.0'], T = 100, Caflux = 100.0, Lflux = 5.0, Gluflux = 10.0, altereds = ['_k1x1.0'], freqs = [0.5, 1.0, 2.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 22.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 200.0, 300.0, 500.0], epsfilename = '', givenlabels = [], givencolors = []):
 if len(blockeds) == 1 and len(altereds) == 1:
  return ltdltpcurve(blockeds[0], T, Caflux, Lflux, Gluflux, altereds[0], freqs)
 else:
  DATAS = []
  print("blockeds="+str(blockeds))
  labels = [] 
  if len(blockeds) > 1 and len(altereds) == 1:
   for iexp in range(0,len(blockeds)):
    DATA = ltdltpcurve(blockeds[iexp], T, Caflux, Lflux, Gluflux, altereds[0], freqs)
    DATAS.append(DATA[:])
    labels.append('blocked-'+str(iexp))
  elif len(blockeds) == 1 and len(altereds) > 1:
   for iexp in range(0,len(altereds)):
    DATA = ltdltpcurve(blockeds[0], T, Caflux, Lflux, Gluflux, altereds[iexp], freqs)
    DATAS.append(DATA[:])
    labels.append('altered-'+str(iexp))
  elif len(blockeds) == len(altereds):
   for iexp in range(0,len(blockeds)):
    DATA = ltdltpcurve(blockeds[iexp], T, Caflux, Lflux, Gluflux, altereds[iexp], freqs)
    DATAS.append(DATA[:])
    labels.append('blocked-'+str(iexp))
  else:
   for iexp in range(0,len(blockeds)):
    for iexp2 in range(0,len(altereds)):
     DATA = ltdltpcurve(blockeds[iexp], T, Caflux, Lflux, Gluflux, altereds[iexp2], freqs)
     DATAS.append(DATA[:])
    labels.append('blocked-'+str(iexp)+'_altered-'+str(iexp2))
  return DATAS 


#Define the simulation attributes
T = 100        #time of stimulation in seconds
Lflux = 5.0    #flux of beta-adrenergic ligands (particles/ms)
Gluflux = 10.0 #flux of glutamate for mGluR activation (particles/ms)

blockeds = ['PP1','PKA','PDE4','PLA2','PKC','NCX','PP1','PKA','PDE4','PLA2','PKC','NCX','PKC','CK','Gqabg','MGluR','Calbin,CalbinC','DAGK','Gi']
blockedTitles = ['PP1','PKA','PDE4','PLA2','PKC','NCX','PP1','PKA','PDE4','PLA2','PKC','NCX','PKC','CaMKII','Gq','mGluR','Calbindin','DAGK','Gi']
toplots = [2,2,2,2,2,2,3,3,3,3,3,3,0,0,0,0,0,0,0]

blockedCoeffTypes = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] # For the proteins, use blockedSets[1] (refers to blockedCoeffs[1] and blockedCoeffs[2], meaning coeffs 0.8 and 1.2; and for combinations, use blockedSets[2] (refers to blockedCoeffs[4] and blockedCoeffs[5], meaning coeffs 10 and 11))
blockedSets = [[0,3],[1,2],[4,5]] 

blockedCoeffs = [0.5,0.8,1.2,1.5,10,11] 

cols = ['#BBBBBB','#666666', '#AA0000', '#FF2222','#FF00FF','#888800']
dimcols = ['#AAAAAA','#999999', '#FF8888','FFDDDD','#FF88FF','#AAAA44']


f,axs = subplots(5,4)
axarr = axs.reshape(prod(axs.shape),).tolist()
axinsets = []
for iax in range(0,len(axarr)):
  for tick in axarr[iax].xaxis.get_major_ticks()+axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(5)
  for axis in ['top','bottom','left','right']:
    axarr[iax].spines[axis].set_linewidth(0.3)
  boxoff(axarr[iax])
for iax in range(0,6):
  axarr[iax].set_position([0.095+0.15*iax,0.8,0.15,0.15])
  axinsets.append(f.add_axes([0.095+0.15*iax+0.04,0.8+0.07,0.01,0.07]))
for iax in range(0,6):
  axarr[6+iax].set_position([0.095+0.15*iax,0.59,0.15,0.15])
  axinsets.append(f.add_axes([0.095+0.15*iax+0.04,0.59+0.07,0.01,0.07]))
for iax in range(0,4):
  axarr[12+iax].set_position([0.095+0.225*iax,0.33,0.225,0.15])
  axinsets.append(f.add_axes([0.095+0.225*iax+0.06,0.33+0.07,0.02,0.1]))
for iax in range(0,4):
  axarr[16+iax].set_position([0.095+0.225*iax,0.1,0.225,0.15])
  axinsets.append(f.add_axes([0.095+0.225*iax+0.06,0.1+0.07,0.02,0.1]))

for iax in range(0,len(axinsets)):
  for tick in axinsets[iax].xaxis.get_major_ticks()+axinsets[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(5)
  for axis in ['top','bottom','left','right']:
    axinsets[iax].spines[axis].set_linewidth(0.3)
  axinsets[iax].set_xticks([])
  boxoff(axinsets[iax])

axarr[19].set_visible(False)
for iax in [0,1,2,3,4,5,6,7,8,9,10,11,19]:
  axinsets[iax].set_visible(False)

for i in [1,2,3,4,5,7,8,9,10,11,13,14,15,17,18,19]:
  axarr[i].set_yticklabels([])
for i in range(0,19):
  points = axarr[i].get_position().get_points()
  f.text(points[0,0]+0.007-0.0*(i<12),points[1,1]-0.024+0.03*(i<12),chr(ord('A')+i),fontsize=8)

### Population-average simulations. Panels A-G.

GluRCoeff = '0.5,0.5,1.5,1.5'
myfreqs = [0.001]+[1.0*i for i in range(1,28)]

Caflux = 100.0

DATAS_0_CONTROL = ltdltpcurves_many(['GluR1,GluR1_memb,GluR2,GluR2_memb,Cax'+GluRCoeff+',1.0'], T, Caflux, Lflux, Gluflux, [''], freqs = myfreqs)
#DATAS_0_CONTROL[0][i][0] contains the post-30min conductance for frequency i for control synapse, and DATAS_0_CONTROL[0][i][1] contains the baseline conductance (should be the same for all frequencies since taken just before the stimulus onset)
minmaxes = []
  
for iblocked in range(0,len(blockeds)):
  blockedCoeffType = blockedCoeffTypes[iblocked]
  blockedSet = blockedSets[blockedCoeffType]
  toplot = toplots[iblocked]
  for iicoeff in range(0,len(blockedSet)):
    iblockedCoeff = blockedSet[iicoeff]
    coeff = blockedCoeffs[iblockedCoeff]
    Nsamespecies = blockeds[iblocked].count(',')+1
    blocked_base = 'GluR1,GluR1_memb,GluR2,GluR2_memb,'+blockeds[iblocked]+'x'+GluRCoeff+','+','.join([str(coeff) for i in range(0,Nsamespecies)])
    DATAS_0 = ltdltpcurves_many([blocked_base], T, Caflux, Lflux, Gluflux, [''], freqs = myfreqs)
    #DATAS_0[0][i][0] contains the post-30min conductance for frequency i for the variant synapse, and DATAS_0[0][i][1] contains the baseline conductance (should be the same for all frequencies since taken just before the stimulus onset)

    axarr[iblocked].plot(myfreqs,[DATAS_0[toplot][i][0] for i in range(0,len(DATAS_0[0]))],'b.-', lw=0.3, ms=1.0, mew=1.0, color=cols[iblockedCoeff], label = '-20%' if iblockedCoeff==1 else '+20%' )
    axinsets[iblocked].bar(-1+2*iicoeff if iblocked < len(blockeds)-1 else iicoeff,DATAS_0[1][0][1],color=cols[iblockedCoeff])
  axarr[iblocked].plot(myfreqs,[DATAS_0_CONTROL[toplot][i][0] for i in range(0,len(DATAS_0_CONTROL[0]))],'b.-', lw=0.3, ms=1.0, mew=1.0, color='#000000', label = 'control')
  axarr[iblocked].set_xlim([0,27])
  #axarr[iblocked].set_ylim([0.6,2.6])

  axinsets[iblocked].bar(0 if iblocked < len(blockeds)-1 else -1,DATAS_0_CONTROL[1][0][1],color='#000000')
  axinsets[iblocked].set_xlim([-1.5,1.5])
  axinsets[iblocked].set_ylim([0,43])

  if iblocked < 6 or iblocked >= 12:
    axarr[iblocked].set_title(blockedTitles[iblocked],fontsize=6.5)

  axarr[0].set_ylabel('[GluR1]$_{\mathrm{memb}}$ (mM)',fontsize=6.5)
  axarr[6].set_ylabel('[GluR2]$_{\mathrm{memb}}$ (mM)',fontsize=6.5)
  axarr[12].set_ylabel('Rel. cond.\n(fold change)',fontsize=6.5)
  axarr[16].set_ylabel('Rel. cond.\n(fold change)',fontsize=6.5)
    
#for iax in [0,4]:
#  axarr[iax].set_ylabel('Rel. cond. (fold change)',fontsize=6.5)
for iax in [6,7,8,9,10,11,15,16,17,18]:
  axarr[iax].set_xlabel('Freq. (Hz)',fontsize=6.5)

  
#handles, labels = axarr[5].get_legend_handles_labels()
#leg = axarr[5].legend([handles[idx] for idx in [0,2,1]],[labels[idx] for idx in [0,2,1]],fontsize=6,frameon=False,bbox_to_anchor=(1.05,0.94))
handles, labels = axarr[6].get_legend_handles_labels()
leg = axarr[6].legend([handles[idx] for idx in [0,2,1]],[labels[idx] for idx in [0,2,1]],fontsize=6,frameon=False,loc=1)
#handles, labels = axarr[0].get_legend_handles_labels()
#leg = axarr[0].legend([handles[idx] for idx in [0,2,1]],[labels[idx] for idx in [0,2,1]],fontsize=6,frameon=False,loc=4)



#f.set_size_inches(8.27,11.69)
f.savefig('fig_ltdltpcurves_mut_all_supp.eps')
