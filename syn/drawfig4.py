from pylab import *
from os.path import exists
import scipy.io
import calcconds
import scipy.stats
from matplotlib.collections import PatchCollection

def drawdiscontinuity(ax,y,yoffset,x=0,xoffset=0.1,lw=2.0,lw2=1.0):
  thisline = ax.plot([x-xoffset,x+xoffset],[y-yoffset,y],'k-',linewidth=lw2)
  thisline[0].set_clip_on(False)
  thisline = ax.plot([x-xoffset,x+xoffset],[y,y+yoffset],'k-',linewidth=lw2)
  thisline[0].set_clip_on(False)
  thisline = ax.plot([x-1.5*xoffset,x+1.5*xoffset],[y-0.75*yoffset,y+0.75*yoffset],'k-',color='#FFFFFF',zorder=100,linewidth=lw)
  thisline[0].set_clip_on(False)

#boxoff: remove extra axes (top and right) from an axis                                                                                                                                                                                                                   
def boxoff(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()


Ntrains = list(range(2,19)) #[2,4,6,8,10,12,14,16,18,20,22,24]
GluRcoeffs = ['2.0,2.0,0.0,0.0', '0.0,0.0,2.0,2.0']

Nstim_MMN = 10
Nstim_VEP = 1

f,axs = subplots(8,2)
axarr = axs.reshape(prod(axs.shape),).tolist()
for iax in range(0,4):
  axarr[iax].set_position([0.08+0.24*iax,0.85,0.17,0.105])
  axarr[iax].set_xlim([0.0,23])
  if iax < 2:
    axarr[iax].set_ylim([0.5,2.2])
  else:
    axarr[iax].set_ylim([0.5,1.06])
axarr[-1].set_visible(False)
    
for iax in range(0,4):
  axarr[4+iax].set_position([0.08+0.23*iax,0.59,0.20,0.19])
for iax in range(0,3):
  #axarr[8+iax].set_position([0.18+0.23*iax,0.59-0.25,0.21,0.20])
  axarr[8+iax].set_position([0.08+0.23*(iax+(iax==2)),0.59-0.25,0.21,0.19])
# axarr[4+iax].set_xlim([0,60])
# axarr[4+iax].set_ylim([0.4,2.7])
for iax in range(0,2):
  axarr[11+iax].set_position([0.08+0.21*iax,0.08,0.19,0.17])
  axarr[13+iax].set_position([0.58+0.21*iax,0.08,0.19,0.17])


for iax in range(0,len(axarr)):
  boxoff(axarr[iax])
  axarr[iax].set_xlim([0,60])
  axarr[iax].set_ylim([0.4,2.4])
  for tick in axarr[iax].xaxis.get_major_ticks()+axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(5)
  for axis in ['top','bottom','left','right']:
    axarr[iax].spines[axis].set_linewidth(0.3)
  axarr[iax].xaxis.set_tick_params(width=0.5)
  axarr[iax].yaxis.set_tick_params(width=0.5)

for iax in [11,12]:
  axarr[iax].set_ylim([1.0,5.7])
for iax in [13,14]:
  axarr[iax].set_ylim([0.25,1.0])
#VEP with NM
styles = ['k-','k--']; dasheds = [(2,0),(1.8,1.8)]
Caflux = 150.0; Lflux = 5.0; Gluflux = 10.0
conds_control = []
times_control = []
for iGluRcoeff in [0,1]:
  filename = 'nrn_tstop32000000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_membx'+GluRcoeffs[iGluRcoeff]+'_k1x1.0_onset24040000.0_n'+str(Nstim_VEP)+'_freq100.0_dur3.0_flux'+str(Caflux)+'_Lflux'+str(Lflux)+'_Gluflux'+str(Gluflux)+'_AChflux'+str(Gluflux)+'_Ntrains1200_trainT500.0.mat'
  if not exists(filename):
    print(filename+' does not exist')
    continue
  print("Loading "+filename)
  conds,times = calcconds.calcconds_nrn(filename)
  axarr[2*iGluRcoeff+1].plot([(x-24040000)/60000 for x in times],[x/conds[0] for x in conds],styles[iGluRcoeff],dashes=dasheds[iGluRcoeff], lw=0.5, label='with NM')
  conds_control.append(conds[:])
  times_control.append(times[:])
#print(str(min(times_control[0])))
#print(str(max(times_control[0])))
#print(str(min(times_control[1])))
#print(str(max(times_control[1])))
#print("i_3600s = "+str([i for i in range(0,len(times_control[0])) if abs(times_control[0][i] - (24040000+3600000)) < 100.0][0]))
#print("t_3600s = "+str(times_control[0][[i for i in range(0,len(times_control[0])) if abs(times_control[0][i] - (24040000+3600000)) < 100.0][0]]))
#qwer

#VEP without NM
Caflux = 150.0; Lflux = 0.0; Gluflux = 0.0
for iGluRcoeff in [0,1]:
  filename = 'nrn_tstop32000000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_membx'+GluRcoeffs[iGluRcoeff]+'_k1x1.0_onset24040000.0_n'+str(Nstim_VEP)+'_freq100.0_dur3.0_flux'+str(Caflux)+'_Lflux'+str(Lflux)+'_Gluflux'+str(Gluflux)+'_AChflux'+str(Gluflux)+'_Ntrains1200_trainT500.0.mat'
  if not exists(filename):
    print(filename+' does not exist')
    continue
  print("Loading "+filename)
  conds,times = calcconds.calcconds_nrn(filename)
  axarr[2*iGluRcoeff].plot([(x-24040000)/60000 for x in times],[x/conds[0] for x in conds],styles[iGluRcoeff],dashes=dasheds[iGluRcoeff],color='#0000FF', lw=0.5, label='without NM')
  print('  VEP amplitude relative to control: '+str(conds_control[iGluRcoeff][361]/conds_control[iGluRcoeff][0])+'=>'+str(conds[361]/conds[0])+' = '+str((conds[361]/conds[0]-1)/(conds_control[iGluRcoeff][361]/conds_control[iGluRcoeff][0]-1)))


#blockedSets = [[0,3],[1,2],[4,7],[5,6],[8,9],[10,11]] #0: +-50%, 1: +-20%, 2: LTP-decreasing and LTD-decreasing, 3: LTP-increasing and LTD-increasing, 4: MMN-decreasing & increasing, 5: VEP-decreasing & increasing
blockedCoeffTypes = [1,1,1,1,1,1,1,1,1,1,1,1,1,2,2]
blockedSets = [[0,3],[1,2],[10,11]] 
blockedCoeffs = [0.5,0.8,1.2,1.5,4,5,6,7,8,9,10,11,12,13] #0: +-50%, 1: +-20%, 2: LTP-decreasing and LTD-decreasing, 3: LTP-increasing and LTD-increasing, 4-5: CM combs from non-imputed-Both, 6-7: CM combs from imputed-GWAS, 8-9: CM combs from imputed-CM, 10-11: CM combs fr#cols = ['#BBBBBB','#666666', '#AA0000', '#FF2222','#660066','#006600','#BB00BB','#00BB00','#AAAA00','#FFFF00','#FF9B00','#000066','#FF9B00','#333388']
#cols = ['#BBBBBB','#666666', '#AA0000', '#FF2222','#660066','#006600','#BB00BB','#00BB00','#AAAA00','#FFFF00','#FF9B00','#000066','#770077','#00FF00']
cols = ['#BBBBBB','#666666', '#AA0000', '#FF2222','#FF00FF','#888800','#FF00FF','#888800','#FF00FF','#888800','#FF00FF','#888800','#FF00FF','#888800','#FF00FF','#888800','#FF00FF','#888800','#FF00FF','#888800']


#VEP with muts
blockeds = ['PP1','PKA','PDE4','PLA2','Ng','NCX','CMcomb','CMcomb']
blockedTitles = ['PP1','PKA','PDE4','PLA2','Neurogranin','NCX','Combinations','Combinations']
blockedCoeffTypes = [1,1,1,1,1,1,2]

Caflux = 150.0; Lflux = 5.0; Gluflux = 10.0
iGluRcoeff = 1
for iblocked in range(0,len(blockedCoeffTypes)):
  blockedCoeffType = blockedCoeffTypes[iblocked]
  blockedSet = blockedSets[blockedCoeffType]
  for iicoeff in range(0,len(blockedSet)):
    iblockedCoeff = blockedSet[iicoeff]
    coeff = blockedCoeffs[iblockedCoeff]
    for iGluRcoeff in [0,1]:
      filename = 'nrn_tstop32000000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,'+blockeds[iblocked]+'x'+GluRcoeffs[iGluRcoeff]+','+str(coeff)+'_k1x1.0_onset24040000.0_n'+str(Nstim_VEP)+'_freq100.0_dur3.0_flux'+str(Caflux)+'_Lflux'+str(Lflux)+'_Gluflux'+str(Gluflux)+'_AChflux'+str(Gluflux)+'_Ntrains1200_trainT500.0.mat'
      if not exists(filename):
        print(filename+' does not exist')
        continue
      print("Loading "+filename)
      conds,times = calcconds.calcconds_nrn(filename)
      addition = ' (GluR1 synapse)' if iGluRcoeff == 0 else ' (GluR2 synapse)'
      if iGluRcoeff == 0 and iblocked != 5 or iGluRcoeff == 1 and iblocked == 5: # Opposite curves (dashed) labeled in axarr[19] than the rest, so they can be printed in an extra legend
        axarr[4+iblocked].plot([(x-24040000)/60000 for x in times],[x/conds[0] for x in conds],styles[iGluRcoeff],dashes=dasheds[iGluRcoeff],color=cols[iblockedCoeff],lw=0.5,label=('-20%' if iblockedCoeff==1 else '+20%')+addition if iblocked < 6 else ('PFC' if iicoeff==0 else 'ACC'))
      else:
        axarr[4+iblocked].plot([(x-24040000)/60000 for x in times],[x/conds[0] for x in conds],styles[iGluRcoeff],dashes=dasheds[iGluRcoeff],color=cols[iblockedCoeff],lw=0.5)
      print('  VEP amplitude relative to control: '+str(conds_control[iGluRcoeff][361]/conds_control[iGluRcoeff][0])+'=>'+str(conds[361]/conds[0])+' = '+str((conds[361]/conds[0]-1)/(conds_control[iGluRcoeff][361]/conds_control[iGluRcoeff][0]-1))+'. times[361]='+str(times[361]))
      #print('    conds_control[iGluRcoeff] = '+str(conds_control[iGluRcoeff][0:380:20])+', conds = '+str(conds[0:380:20]))
  for iGluRcoeff in [0,1]:
    addition = (' (GluR1 synapse)' if iGluRcoeff == 0 else ' (GluR2 synapse)') if iblocked < 6 else ''
    if iGluRcoeff == 0 and iblocked != 5 or iGluRcoeff == 1 and iblocked == 5:
      axarr[4+iblocked].plot([(x-24040000)/60000 for x in times_control[iGluRcoeff]],[x/conds_control[iGluRcoeff][0] for x in conds_control[iGluRcoeff]],styles[iGluRcoeff],dashes=dasheds[iGluRcoeff],lw=0.5,label='control'+addition)
    else:
      axarr[4+iblocked].plot([(x-24040000)/60000 for x in times_control[iGluRcoeff]],[x/conds_control[iGluRcoeff][0] for x in conds_control[iGluRcoeff]],styles[iGluRcoeff],dashes=dasheds[iGluRcoeff],lw=0.5)
  #axarr[4+iblocked].set_title(blockedTitles[iblocked],fontsize=5)
  axarr[4+iblocked].text(5,2.4+0.12*(iblocked==9),blockedTitles[iblocked],fontsize=6,ha='left')

for iax in range(0,len(axarr)-1):
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.02, pos.y1 + 0.01, chr(ord('A')+iax), fontsize=8)

#for iax in [4,5,6,7,8]:
#  axarr[iax].set_xticklabels([])
for iax in [0,1,2,3,4,5,6,7,8,9,10]:
  axarr[iax].set_xlabel('time (min)',fontsize=6.5,labelpad=1)

axarr[0].set_title('GluR1 synapse', fontsize=6.5)
axarr[1].set_title('GluR1 synapse', fontsize=6.5)
axarr[2].set_title('GluR2 synapse', fontsize=6.5)
axarr[3].set_title('GluR2 synapse', fontsize=6.5)
#f.text(0.02,0.16,'VEP, with NM', fontsize=6, rotation=90,ha='right',va='center')

axarr[0].legend(fontsize=5.6,frameon=False)
axarr[1].legend(fontsize=5.6,frameon=False)
axarr[2].legend(fontsize=5.6,frameon=False)
axarr[3].legend(fontsize=5.6,frameon=False)

handles, labels = axarr[6].get_legend_handles_labels()
leg = axarr[6].legend([handles[idx] for idx in [0,2,1]],[labels[idx] for idx in [0,2,1]],fontsize=5.6,frameon=False,bbox_to_anchor=(1.01,-0.27))
handles, labels = axarr[9].get_legend_handles_labels()
leg = axarr[9].legend([handles[idx] for idx in [0,2,1]],[labels[idx] for idx in [0,2,1]],fontsize=5.6,frameon=False,bbox_to_anchor=(1.03,0.29))
axarr[10].legend(fontsize=5.6,frameon=False,bbox_to_anchor=(0.35,0.35))

f.savefig('figVEP.eps')

commonattrs = 'isSign0_isImputed1_Both'
pvalthresh = 0.05


def drawdiscontinuity(ax,y,yoffset,x=0,xoffset=0.1,lw=2.0,lw2=1.0):
  thisline = ax.plot([x-xoffset,x+xoffset],[y-yoffset,y],'k-',linewidth=lw2)
  thisline[0].set_clip_on(False)
  thisline = ax.plot([x-xoffset,x+xoffset],[y,y+yoffset],'k-',linewidth=lw2)
  thisline[0].set_clip_on(False)
  thisline = ax.plot([x-1.5*xoffset,x+1.5*xoffset],[y-0.75*yoffset,y+0.75*yoffset],'k-',color='#FFFFFF',zorder=100,linewidth=lw)
  thisline[0].set_clip_on(False)

def mybar(ax,x,y,facecolor=[],linewidth=0.3,w=0.4):
  if len(y) == 0:
    return [nan,nan]
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

def mysigndiff(ax,x0,x1,y0,y1):
  #maxval = max(mean(y0)+std(y0),mean(y1)+std(y1))
  maxval = max(max(y0),max(y1))
  ax.plot([x0,x0,x1,x1],[maxval*1.017,maxval*1.034,maxval*1.034,maxval*1.017],'k-',lw=0.3)
  ax.plot(0.5*(x0+x1),maxval*1.06,'k*',ms=1.8,mew=0.3,lw=0.3)
def mysigndiff2(ax,x0,x1,y0,y1):
  #maxval = max(mean(y0)+std(y0),mean(y1)+std(y1))
  maxval = max(max(y0),max(y1))
  ax.plot(0.5*(x0+x1),maxval*1.09,'r*',ms=1.8,mew=0.3,lw=0.3)
  
Ntrain = 18
GluRcoeffs = ['2.0,2.0,0.0,0.0', '0.0,0.0,2.0,2.0']

Nstim_MMN = 10
Nstim_VEP = 1

areas = ['PFC','ACC']
for iax in range(0,4):
  axarr[11+iax].set_title(areas[iax%2]+', GluR'+str(1+int(iax/2)),fontsize=6.5)
axarr[0].set_ylabel('Relative syn. strength\n(fold change)',fontsize=6.5)
axarr[4].set_ylabel('Relative syn. strength\n(fold change)',fontsize=6.5)
axarr[8].set_ylabel('Relative syn. strength\n(fold change)',fontsize=6.5)
axarr[11].set_ylabel('Relative syn. strength\n(fold change)',fontsize=6.5)
axarr[13].set_ylabel('Relative syn. strength\n(fold change)',fontsize=6.5)

cols = ['#666666', '#FF00FF','#888800']
dimcols = ['#AAAAAA','#FF88FF','#AAAA44']

#MMN (LTP+LTD) with NM
Caflux = 150.0; Lflux = 5.0; Gluflux = 10.0
areas = ['PFC','ACC']
for iarea in [0,1]:
  area = areas[iarea]
  attrs = area + '_'+commonattrs
  groups = ['HC','SCZ']
  pvals_all = []
  datas = [[],[]]
  pvals = []
  #From extract_commonmind_ACC_data_sklearn_allgenes_params_samples.py:
  if area == 'ACC':
    isubjs_HC = [1, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 73, 75, 76, 77, 78, 79, 80, 83, 84, 86, 90, 91, 92, 95, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 124, 126, 127, 129, 130, 131, 132, 133, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 152, 154, 159, 161, 163, 164, 167, 169, 170, 172, 174, 175, 176, 177, 182, 199, 208, 235, 240, 243, 249, 253, 256, 257, 260, 261, 262, 263, 265, 266, 267, 268, 271, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 329, 330, 331, 332, 343, 349, 350, 352, 353, 355, 360, 366, 369, 370, 371, 372, 374, 375, 376, 378, 381, 382, 383, 385, 388, 389, 392, 393, 395, 396, 397, 398, 401, 402, 409, 410, 414, 416, 418, 421, 422, 426, 427, 431, 432, 433, 434, 435, 437, 438, 442, 443, 444, 446, 447, 448, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480]
    isubjs_SCZ = [0, 2, 3, 12, 17, 28, 31, 58, 69, 70, 71, 72, 74, 81, 82, 85, 87, 88, 89, 93, 94, 96, 97, 98, 99, 100, 116, 117, 118, 119, 120, 121, 122, 123, 125, 134, 151, 155, 156, 157, 158, 160, 162, 165, 166, 168, 171, 173, 178, 179, 180, 181, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 200, 201, 202, 203, 204, 205, 206, 207, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 236, 237, 238, 239, 241, 242, 244, 245, 246, 247, 248, 250, 251, 252, 254, 255, 258, 259, 264, 269, 270, 272, 273, 274, 275, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 344, 345, 346, 347, 348, 351, 354, 356, 357, 358, 359, 361, 362, 363, 364, 365, 367, 368, 373, 377, 379, 380, 384, 386, 387, 390, 391, 394, 399, 400, 403, 404, 405, 406, 407, 411, 412, 413, 415, 417, 419, 420, 423, 424, 425, 428, 429, 430, 436, 439, 440, 441, 445, 449, 450, 451, 452]
  elif area == 'PFC':
    #From extract_commonmind_data_sklearn_allgenes_params_samples.py:
    isubjs_HC = [0, 1, 2, 3, 5, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 19, 20, 22, 23, 25, 27, 28, 31, 32, 33, 34, 35, 36, 37, 39, 41, 42, 43, 47, 51, 52, 53, 54, 63, 64, 66, 69, 73, 80, 81, 82, 83, 84, 87, 88, 89, 90, 92, 93, 94, 96, 97, 100, 105, 107, 111, 112, 118, 119, 120, 125, 127, 128, 135, 139, 140, 142, 143, 145, 146, 152, 154, 156, 165, 168, 176, 181, 182, 187, 192, 194, 201, 203, 204, 205, 207, 210, 211, 217, 220, 222, 223, 225, 226, 232, 235, 236, 240, 241, 242, 243, 244, 247, 248, 249, 251, 252, 253, 254, 257, 258, 262, 263, 265, 266, 267, 268, 269, 270, 271, 272, 283, 290, 291, 293, 294, 295, 296, 297, 298, 299, 300, 301, 305, 306, 308, 313, 314, 315, 317, 319, 320, 321, 323, 324, 325, 326, 329, 330, 332, 333, 335, 336, 339, 340, 341, 342, 343, 344, 345, 347, 351, 352, 354, 355, 356, 359, 361, 362, 368, 371, 373, 374, 375, 376, 377, 380, 381, 382, 383, 384, 386, 387, 388, 389, 390, 393, 396, 397, 398, 399, 400, 402, 403, 405, 406, 407, 408, 415, 416, 417, 420, 421, 422, 423, 424, 425, 426, 427, 428]
    isubjs_SCZ = [4, 6, 15, 17, 21, 24, 26, 29, 30, 38, 40, 44, 45, 46, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, 62, 65, 67, 68, 70, 71, 72, 74, 75, 76, 77, 78, 79, 85, 86, 91, 95, 98, 99, 101, 102, 103, 104, 106, 108, 109, 110, 113, 114, 115, 116, 117, 121, 122, 123, 124, 126, 129, 130, 131, 132, 133, 134, 136, 137, 138, 141, 144, 147, 148, 149, 150, 151, 153, 155, 157, 158, 159, 160, 161, 162, 163, 164, 166, 167, 169, 170, 171, 172, 173, 174, 175, 177, 178, 179, 180, 183, 184, 185, 186, 188, 189, 190, 191, 193, 195, 196, 197, 198, 199, 200, 202, 206, 208, 209, 212, 213, 214, 215, 216, 218, 219, 221, 224, 227, 228, 229, 230, 231, 233, 234, 237, 238, 239, 245, 246, 250, 255, 256, 259, 260, 261, 264, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 285, 286, 287, 288, 289, 292, 302, 303, 304, 307, 309, 310, 311, 312, 316, 318, 322, 327, 328, 334, 337, 338, 346, 348, 349, 350, 353, 357, 358, 360, 363, 364, 365, 366, 367, 369, 370, 372, 378, 379, 385, 392, 394, 395, 401, 404, 409, 410, 411, 412, 413, 414, 418, 419]
  else:
    error('Unknown area')
  isubjs_HC = [0,1,2,3,4]; isubjs_SCZ = [5,6,7,8,9]
    
  #VEP with NM
  styles = ['k-','k--']; dasheds = [(2,0),(1.8,1.8)]
  Caflux = 150.0; Lflux = 5.0; Gluflux = 10.0
  pvals_all = []
  for iglur in [0,1]:
    pvals_thisglur = []
    datas_thisglur = [[],[]]
    for isamp in range(0,481):
        if isamp in isubjs_HC:
          igroup = 0
        elif isamp in isubjs_SCZ:
          igroup = 1
        else:
          continue
        filename = 'nrn_tstop32000000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,CMcomb_samps_nonIC_'+attrs+'x'+GluRcoeffs[iglur]+','+str(isamp)+'_k1x1.0_onset24040000.0_n'+str(Nstim_VEP)+'_freq100.0_dur3.0_flux'+str(Caflux)+'_Lflux'+str(Lflux)+'_Gluflux'+str(Gluflux)+'_AChflux'+str(Gluflux)+'_Ntrains1200_trainT500.0.mat'
        if not exists(filename):
          print(filename+' does not exist')
          continue
        if isamp % 100 == 99:
          print("Loading "+filename+" ("+str(isamp)+" done)")
        conds,times = calcconds.calcconds_nrn(filename)
        datas_thisglur[igroup].append([x/conds[0] for x in conds])
        #print(str(len(datas_thisglur[0]))+' '+str(len(datas_thisglur[1])))
    for it in range(0,11*60): #len(times)):
        pval = scipy.stats.ranksums(array([datas_thisglur[0][isamp][it] for isamp in range(0,len(datas_thisglur[0]))]),array([datas_thisglur[1][isamp][it] for isamp in range(1,len(datas_thisglur[1]))]))[1]
        #print('iax='+str(iax)+', it='+str(it)+', pval='+str(pval))
        pvals.append(pval)
        if it%60 == 1:
          mybar(axarr[11+2*iglur+iarea],2*int(it/60),[datas_thisglur[0][isamp][it] for isamp in range(0,len(datas_thisglur[0]))],facecolor='#666666',w=0.375)
          mybar(axarr[11+2*iglur+iarea],2*int(it/60)+0.75,[datas_thisglur[1][isamp][it] for isamp in range(0,len(datas_thisglur[1]))],facecolor=cols[1+iarea],w=0.375)
          print("it = "+str(it)+", pval = "+str(pval))
          if pval < pvalthresh and it > 1:
            mysigndiff(axarr[11+2*iglur+iarea],2*int(it/60),2*int(it/60)+0.75,[datas_thisglur[0][isamp][it] for isamp in range(0,len(datas_thisglur[0]))],[datas_thisglur[1][isamp][it] for isamp in range(0,len(datas_thisglur[1]))])
            if pval < pvalthresh/11:
              mysigndiff2(axarr[11+2*iglur+iarea],2*int(it/60),2*int(it/60)+0.75,[datas_thisglur[0][isamp][it] for isamp in range(0,len(datas_thisglur[0]))],[datas_thisglur[1][isamp][it] for isamp in range(0,len(datas_thisglur[1]))])
        
    print('iglur='+str(iglur)+', iarea='+str(iarea)+', min pval='+str(min(pvals)))
    pvals_all.append(pvals[:])

for iax in [0,1]:
  for iglur in [0,1]:
    axarr[11+2*iglur+iax].set_xlabel('time (min)',fontsize=6.5)
    axarr[11+2*iglur+iax].set_xticks([2*i+0.375 for i in range(0,11)])
    axarr[11+2*iglur+iax].set_xticklabels([str(6*i) for i in range(0,11)])
    axarr[11+2*iglur+iax].set_xlim([0,21.5])
    if iax == 1:
      axarr[11+2*iglur+iax].set_yticklabels([])

#axarr[13].legend(fontsize=5.6,frameon=False)
#axarr[14].legend(fontsize=5.6,frameon=False)

polygon = Polygon(array([[12,12,13,13],[0.39,0.44,0.44,0.39]]).T, True)
p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
p.set_facecolor(cols[0])
p.set_edgecolor('None')
axarr[13].add_collection(p)

polygon = Polygon(array([[12,12,13,13],[0.3,0.35,0.35,0.3]]).T, True)
p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
p.set_facecolor(cols[1])
p.set_edgecolor('None')
axarr[13].add_collection(p)

polygon = Polygon(array([[12,12,13,13],[0.39,0.44,0.44,0.39]]).T, True)
p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
p.set_facecolor(cols[0])
p.set_edgecolor('None')
axarr[14].add_collection(p)

polygon = Polygon(array([[12,12,13,13],[0.3,0.35,0.35,0.3]]).T, True)
p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
p.set_facecolor(cols[2])
p.set_edgecolor('None')
axarr[14].add_collection(p)

axarr[13].text(14,0.41,'Control',fontsize=6,ha='left',va='center')
axarr[14].text(14,0.41,'Control',fontsize=6,ha='left',va='center')
axarr[13].text(14,0.32,'SCZ (PFC)',fontsize=6,ha='left',va='center')
axarr[14].text(14,0.32,'SCZ (ACC)',fontsize=6,ha='left',va='center')

f.set_size_inches([6.4,6.4])
f.savefig('figVEP.eps')
