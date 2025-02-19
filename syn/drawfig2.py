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
def ltdltpcurve(blocked='GluR1,GluR1_memb,GluR2,GluR2_memb,Cax0.5,0.5,1.5,1.5,1.0', T = 100, Caflux = 150.0, Lflux = 5.0, Gluflux = 10.0, altered = '_k1x1.0', freqs = [0.5, 1.0, 2.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 22.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 200.0, 300.0, 500.0]):

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
def ltdltpcurves_many(blockeds=['GluR1,GluR1_memb,GluR2,GluR2_memb,Cax0.5,0.5,1.5,1.5,1.0'], T = 100, Caflux = 150.0, Lflux = 5.0, Gluflux = 10.0, altereds = ['_k1x1.0'], freqs = [0.5, 1.0, 2.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 22.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 200.0, 300.0, 500.0], epsfilename = '', givenlabels = [], givencolors = []):
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

blockeds = ['PP1','PKA','PDE4','PLA2','Ng','NCX','CMcomb']
blockedTitles = ['PP1','PKA','PDE4','PLA2','Neurogranin','NCX','Combinations']

blockedCoeffTypes = [1,1,1,1,1,1,2] # For the proteins, use blockedSets[1] (refers to blockedCoeffs[1] and blockedCoeffs[2], meaning coeffs 0.8 and 1.2; and for combinations, use blockedSets[2] (refers to blockedCoeffs[4] and blockedCoeffs[5], meaning coeffs 10 and 11))
blockedSets = [[0,3],[1,2],[4,5]] 

blockedCoeffs = [0.5,0.8,1.2,1.5,10,11] 

cols = ['#BBBBBB','#666666', '#AA0000', '#FF2222','#FF00FF','#888800']
dimcols = ['#AAAAAA','#999999', '#FF8888','FFDDDD','#FF88FF','#AAAA44']


f,axs = subplots(3,5)
axinsets = []
for iax in range(0,5):
  for iay in range(0,3):
    for tick in axs[iay,iax].xaxis.get_major_ticks()+axs[iay,iax].yaxis.get_major_ticks():
      tick.label.set_fontsize(5)
    for axis in ['top','bottom','left','right']:
      axs[iay,iax].spines[axis].set_linewidth(0.3)
    boxoff(axs[iay,iax])

for iay in range(0,3):
  for iax in range(0,5):
    axinsets.append(f.add_axes([0.11+0.18*iax,0.79-0.3*iay,0.03,0.1]))
    for tick in axinsets[-1].xaxis.get_major_ticks()+axinsets[-1].yaxis.get_major_ticks():
      tick.label.set_fontsize(5)
    for axis in ['top','bottom','left','right']:
      axinsets[-1].spines[axis].set_linewidth(0.3)
    axinsets[-1].set_xticks([])
    boxoff(axinsets[-1])

axarr = axs.reshape(prod(axs.shape),).tolist()
for iax in range(0,4):
  axarr[iax].set_position([0.08+0.23*iax,0.74,0.21,0.20])
  axinsets[iax].set_position([0.125+0.23*iax,0.85,0.03,0.08])
for iax in range(0,3):
  axarr[4+iax].set_position([0.18+0.23*iax,0.74-0.31,0.21,0.20])
  axinsets[4+iax].set_position([0.225+0.23*iax,0.85-0.31,0.03,0.08])
for iax in range(0,2):
  axarr[7+iax].set_position([0.08+0.21*iax,0.74-0.31*2,0.19,0.20])
  axarr[9+iax].set_position([0.58+0.21*iax,0.74-0.31*2,0.19,0.20])



for i in range(7,15):
  axinsets[i].set_visible(False)
for i in range(11,15):
  axarr[i].set_visible(False)
for i in [1,2,3,5,6,8,10]:
  axarr[i].set_yticklabels([])
for i in range(0,11):
  points = axarr[i+(i==13)].get_position().get_points()
  f.text(points[0,0]+0.007-0.03*(i>=9),points[1,1]-0.024,chr(ord('A')+i),fontsize=8)

### Population-average simulations. Panels A-G.

GluRCoeffs = ['0.5,0.5,1.5,1.5']
myfreqs = [0.001]+[1.0*i for i in range(1,21)]

Caflux = 150.0
for iglur in range(0,len(GluRCoeffs)):
  DATAS_0_CONTROL = ltdltpcurves_many(['GluR1,GluR1_memb,GluR2,GluR2_memb,Cax'+GluRCoeffs[iglur]+',1.0'], T, Caflux, Lflux, Gluflux, [''], freqs = myfreqs)
  #DATAS_0_CONTROL[0][i][0] contains the post-30min conductance for frequency i for control synapse, and DATAS_0_CONTROL[0][i][1] contains the baseline conductance (should be the same for all frequencies since taken just before the stimulus onset)
  minmaxes = []
  
  for iblocked in range(0,len(blockeds)):
    blockedCoeffType = blockedCoeffTypes[iblocked]
    blockedSet = blockedSets[blockedCoeffType]
    for iicoeff in range(0,len(blockedSet)):
      iblockedCoeff = blockedSet[iicoeff]
      coeff = blockedCoeffs[iblockedCoeff]
      Nsamespecies = blockeds[iblocked].count(',')+1
      blocked_base = 'GluR1,GluR1_memb,GluR2,GluR2_memb,'+blockeds[iblocked]+'x'+GluRCoeffs[iglur]+','+','.join([str(coeff) for i in range(0,Nsamespecies)])
      DATAS_0 = ltdltpcurves_many([blocked_base], T, Caflux, Lflux, Gluflux, [''], freqs = myfreqs)
      #DATAS_0[0][i][0] contains the post-30min conductance for frequency i for the variant synapse, and DATAS_0[0][i][1] contains the baseline conductance (should be the same for all frequencies since taken just before the stimulus onset)

      axarr[iblocked].plot(myfreqs,[DATAS_0[0][i][0] for i in range(0,len(DATAS_0[0]))],'b.-', lw=0.3, ms=1.0, mew=1.0, color=cols[iblockedCoeff], label = ('-20%' if iblockedCoeff==1 else '+20%') if iblocked < 6 else ('PFC' if iicoeff==0 else 'ACC'))
      axinsets[iblocked].bar(-1+2*iicoeff if iblocked < len(blockeds)-1 else iicoeff,DATAS_0[1][0][1],color=cols[iblockedCoeff])
      if iblocked == len(blockeds)-1:
        print(('PFC' if iicoeff==0 else 'ACC')+': max LTP '+str(max([DATAS_0[0][i][0] for i in range(0,len(DATAS_0[0]))])/max([DATAS_0_CONTROL[0][i][0] for i in range(0,len(DATAS_0_CONTROL[0]))]))+' times CTRL - this means '+
              str((max([DATAS_0[0][i][0] for i in range(0,len(DATAS_0[0]))])-1)/(max([DATAS_0_CONTROL[0][i][0] for i in range(0,len(DATAS_0_CONTROL[0]))])-1)*100-100)+'% change in LTP amplitude')

    axarr[iblocked].plot(myfreqs,[DATAS_0_CONTROL[0][i][0] for i in range(0,len(DATAS_0_CONTROL[0]))],'b.-', lw=0.3, ms=1.0, mew=1.0, color='#000000', label = 'Control')
    axarr[iblocked].set_xlim([0,20])
    axarr[iblocked].set_ylim([0.6,2.6])

    axinsets[iblocked].bar(0 if iblocked < len(blockeds)-1 else -1,DATAS_0_CONTROL[1][0][1],color='#000000')
    axinsets[iblocked].set_xlim([-1.5,1.5])
    axinsets[iblocked].set_ylim([0,43])

    axarr[iblocked].set_title(blockedTitles[iblocked],fontsize=6.5)

for iax in [0,4]:
  axarr[iax].set_ylabel('Rel. cond. (fold change)',fontsize=6.5)
for iax in [0,1,2,3,4,5,6,7]:
  axarr[iax].set_xlabel('Freq. (Hz)',fontsize=6.5)

  
handles, labels = axarr[5].get_legend_handles_labels()
leg = axarr[5].legend([handles[idx] for idx in [0,2,1]],[labels[idx] for idx in [0,2,1]],fontsize=6,frameon=False,bbox_to_anchor=(1.05,0.94))
handles, labels = axarr[6].get_legend_handles_labels()
leg = axarr[6].legend([handles[idx] for idx in [0,2,1]],[labels[idx] for idx in [0,2,1]],fontsize=6,frameon=False,loc=4)
handles, labels = axarr[0].get_legend_handles_labels()
leg = axarr[0].legend([handles[idx] for idx in [0,2,1]],[labels[idx] for idx in [0,2,1]],fontsize=6,frameon=False,loc=4)


### Sample-wise simulations. Panels H-K.
freqs = [12.0,14.0,16.0,18.0,20.0]
#freqs = [20.0,24.0,28.0]

freq_to_plot = 16.0
areas = ['PFC','ACC']
attrcommon = 'isSign0_isImputed1_Both'

for iax in range(0,2):
  axarr[0+7+iax].set_xlabel('time (min)',fontsize=6.5)
  axarr[2+7+iax].set_xticks([])
  axarr[0+7+iax].set_title(areas[iax]+", "+str(freq_to_plot)+" Hz", fontsize=6.5)
  axarr[2+7+iax].set_title(areas[iax], fontsize=6.5)
axarr[0+7].set_ylabel('Relative synaptic\nstrength (fold change)',fontsize=6.5)
axarr[2+7].set_ylabel('Relative synaptic strength\nat 30 min (fold change)',fontsize=6.5)

for iax in range(0,4):
  for tick in axarr[7+iax].xaxis.get_major_ticks() + axarr[7+iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(5)

def mybar(ax,x,y,facecolor=[],linewidth=0.3):
  if type(facecolor) is not list or len(facecolor) > 0:
    a1 = ax.bar(x,mean(y),facecolor=facecolor)
  else:
    a1 = ax.bar(x,mean(y))
  a2 = ax.plot([x,x],[mean(y)-std(y),mean(y)+std(y)],'k-',lw=linewidth)
  return [a1,a2]
def mybar2(ax,x,y,facecolor=[],linewidth=0.3,w=0.4):
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
  #isubjs_HC = [0,1,2,3,4]; isubjs_SCZ = [5,6,7,8,9]  #Uncomment this to skip plotting the sample-wise data and replace it with small fake samples

  tpost = 30
  condposts_thisarea = []
  condposts_abs_thisarea = []
  condtots_HC = []
  condtots_SCZ = []

  for ifreq in range(0,len(freqs)):
    freq = freqs[ifreq]
    condposts = []
    condposts_abs = []
    for isamp in range(0,481):
      filename =  'nrn_tstop27000000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,CMcomb_samps_nonIC_'+attr+'x0.5,0.5,1.5,1.5,'+str(isamp)+'_onset24040000.0_n'+str(int(100*freq))+'_freq'+str(freq)+'_dur3.0_flux150.0_Lflux5.0_Gluflux10.0_AChflux10.0.mat'
      if exists(filename):
        print('Loading '+filename)
        A = scipy.io.loadmat(filename)
        conds, times = calcconds.calcconds_nrn(filename)
        minabs = min(abs(times-(24040000+tpost*60*1000)))
        its = [it for it in range(0,len(times)) if times[it]-(24040000+tpost*60*1000) == minabs][0]
        condposts.append(conds[its]/conds[0])
        condposts_abs.append(conds[its])
        if freq == freq_to_plot and isamp in isubjs_HC:
          condtots_HC.append(conds[:])
        if freq == freq_to_plot and isamp in isubjs_SCZ:
          condtots_SCZ.append(conds[:])
      else:
        print(filename+' does not exist')
    if len(condposts) == 0:
      continue
    condposts_thisarea.append(condposts[:])
    condposts_abs_thisarea.append(condposts_abs[:])
    pval = scipy.stats.ranksums(array([condposts[i] for i in isubjs_HC]),array([condposts[i] for i in isubjs_SCZ]))[1]
    print(area+',freq='+str(freq)+',av. condposts='+str(mean(array([condposts[i] for i in isubjs_HC])))+' vs '+str(mean(array([condposts[i] for i in isubjs_SCZ])))+', pval='+str(pval))

    mybar2(axarr[2+7+iarea],3*ifreq,[condposts[i] for i in isubjs_HC],facecolor='#666666')
    mybar2(axarr[2+7+iarea],3*ifreq+1,[condposts[i] for i in isubjs_SCZ],facecolor=cols[4+iarea])
    if pval < 0.05:
      #y1 = 1.03*max(mean([condposts[i] for i in isubjs_HC])+std([condposts[i] for i in isubjs_HC]),mean([condposts[i] for i in isubjs_SCZ])+std([condposts[i] for i in isubjs_SCZ]))
      y1 = 1.03*max(max([condposts[i] for i in isubjs_HC]),max([condposts[i] for i in isubjs_SCZ]))
      axarr[2+7+iarea].plot([3*ifreq,3*ifreq,3*ifreq+1,3*ifreq+1],[y1,y1*1.03,y1*1.03,y1],'k-',lw=0.2)
      axarr[2+7+iarea].plot(3*ifreq+0.5,y1*1.08,'k*',mew=0.05,ms=3.0)
    #axarr[2+7+iarea].text(3*ifreq+0.5,3.0+0.2*iarea,str(freq)+" Hz",ha='center',va='bottom',rotation=90,fontsize=6)
    axarr[2+7+iarea].set_ylim([1.25,3.2])
    axarr[2+7+iarea].set_yticks([1.5,2.0,2.5])
  axarr[2+7+iarea].set_xticks([3*ifreq+0.5 for ifreq in range(0,len(freqs))])
  axarr[2+7+iarea].set_xticklabels([str(int(x)) for x in freqs])
  axarr[2+7+iarea].set_xlabel('Freq. (Hz)',fontsize=6.5)
  
  #condposts_all.append(condposts_thisarea[:])
  #condposts_abs_all.append(condposts_abs_thisarea[:])
  #for ifreq in range(0,len(freqs)):
  #  pval = scipy.stats.ranksums(array([condposts_all[0][ifreq]),array())[1]

  for isamp in range(0,len(condtots_HC)):
    axarr[0+7+iarea].plot((times-24040000)/60/1000,condtots_HC[isamp]/condtots_HC[isamp][0],color=dimcols[0],lw=0.2)
  for isamp in range(0,len(condtots_SCZ)):
    axarr[0+7+iarea].plot((times-24040000)/60/1000,condtots_SCZ[isamp]/condtots_SCZ[isamp][0],color=dimcols[4+iarea],lw=0.2)
  if len(condtots_HC) > 0:
    axarr[0+7+iarea].plot((times-24040000)/60/1000,[mean([condtots_HC[isamp][i]/condtots_HC[isamp][0] for isamp in range(0,len(condtots_HC))]) for i in range(0,len(condtots_HC[0]))],color='#000000',label='Control')
  if len(condtots_SCZ) > 0:
    axarr[0+7+iarea].plot((times-24040000)/60/1000,[mean([condtots_SCZ[isamp][i]/condtots_SCZ[isamp][0] for isamp in range(0,len(condtots_SCZ))]) for i in range(0,len(condtots_SCZ[0]))],color=cols[4+iarea],label='SCZ')
  axarr[0+7+iarea].legend(fontsize=6,frameon=False)

#f.set_size_inches(8.27,11.69)
f.savefig('fig2.eps')
