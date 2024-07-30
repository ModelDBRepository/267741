#cp drawfig_ltdltpcurves_muts_all_ramakercomb.py drawfig_ltdltpcurves_muts_all.py #7.11.2022
from pylab import *
import scipy.io
from os.path import exists
import os
import mytools
import calcconds
from matplotlib.collections import PatchCollection

def boxoff(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def ltdltpcurve(blocked='GluR1,GluR1_memb,GluR2,GluR2_memb,Cax0.5,0.5,1.5,1.5,1.0', T = 100, Caflux = 150.0, Lflux = 5.0, Gluflux = 10.0, altered = '_k1x1.0', freqs = [0.5, 1.0, 2.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 22.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 200.0, 300.0, 500.0],dodraw = True):

 if dodraw:
  close("all")
  f,axarr = subplots(3,4)
  for iay in range(0,3):
   for iax in range(0,4):
     axarr[iay,iax].set_position([0.08+0.25*iax,0.77-0.27*iay,0.17,0.21])
 Nstims = [int(T*x) for x in freqs]

 tposts = [20, 30, 40, 0] #min after stimulus onset

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
    #print("max(conds) = "+str(max(conds))+", min(conds) = "+str(min(conds)))
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
  
def ltdltpcurves_many(blockeds=['GluR1,GluR1_memb,GluR2,GluR2_memb,Cax0.5,0.5,1.5,1.5,1.0'], T = 100, Caflux = 150.0, Lflux = 5.0, Gluflux = 10.0, altereds = ['_k1x1.0'], freqs = [0.5, 1.0, 2.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 22.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 200.0, 300.0, 500.0],dodraw = True, epsfilename = '', givenlabels = [], givencolors = []):
 tposts = [20, 30, 40] #min after stimulus onset
 if len(blockeds) == 1 and len(altereds) == 1:
  return ltdltpcurve(blockeds[0], T, Caflux, Lflux, Gluflux, altereds[0], freqs, dodraw)
 else:
  DATAS = []
  print("blockeds="+str(blockeds))
  labels = [] 
  if len(blockeds) > 1 and len(altereds) == 1:
   for iexp in range(0,len(blockeds)):
    DATA = ltdltpcurve(blockeds[iexp], T, Caflux, Lflux, Gluflux, altereds[0], freqs, False)
    DATAS.append(DATA[:])
    labels.append('blocked-'+str(iexp))
  elif len(blockeds) == 1 and len(altereds) > 1:
   for iexp in range(0,len(altereds)):
    DATA = ltdltpcurve(blockeds[0], T, Caflux, Lflux, Gluflux, altereds[iexp], freqs, False)
    DATAS.append(DATA[:])
    labels.append('altered-'+str(iexp))
  elif len(blockeds) == len(altereds):
   for iexp in range(0,len(blockeds)):
    DATA = ltdltpcurve(blockeds[iexp], T, Caflux, Lflux, Gluflux, altereds[iexp], freqs, False)
    DATAS.append(DATA[:])
    labels.append('blocked-'+str(iexp))
  else:
   for iexp in range(0,len(blockeds)):
    for iexp2 in range(0,len(altereds)):
     DATA = ltdltpcurve(blockeds[iexp], T, Caflux, Lflux, Gluflux, altereds[iexp2], freqs, False)
     DATAS.append(DATA[:])
    labels.append('blocked-'+str(iexp)+'_altered-'+str(iexp2))
  return DATAS 




T = 100
Lflux = 5.0
Gluflux = 10.0

cols = ['#000000','#770077','#009900']


f,axs = subplots(4,6)
axinsets = []
for iax in range(0,6):
  for iay in range(0,4):
    axs[iay,iax].set_position([0.06+0.155*iax,0.065+0.24*(3-iay),0.155,0.205])
    for tick in axs[iay,iax].xaxis.get_major_ticks()+axs[iay,iax].yaxis.get_major_ticks():
      tick.label.set_fontsize(4)
    for axis in ['top','bottom','left','right']:
      axs[iay,iax].spines[axis].set_linewidth(0.3)
    axs[iay,0].set_ylabel('Syn. cond. (fold change)',fontsize=5)
    if (iax == 2 or iax == 3) and (iay == 2 or iay == 3):
      polygon = Polygon(array([[0,40,40,0],[0,0,4,4]]).T, True)
      p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
      p.set_facecolor('#EEEEEE')
      p.set_edgecolor('none')
      axs[iay,iax].add_collection(p)
    if (iax == 4 or iax == 5) and (iay == 1):
      polygon = Polygon(array([[0,40,40,0],[0,0,4,4]]).T, True)
      p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
      p.set_facecolor('#FFE8E8')
      p.set_edgecolor('none')
      axs[iay,iax].add_collection(p)
  axs[3,iax].set_xlabel('$f$ (Hz)',fontsize=5)

for iay in range(0,4):
  for iax in range(0,6):
    axinsets.append(f.add_axes([0.10+0.155*iax,0.117+0.24*(3-iay),0.02,0.09]))         
    for tick in axinsets[-1].xaxis.get_major_ticks()+axinsets[-1].yaxis.get_major_ticks():                                                                                                                                                                           
      tick.label.set_fontsize(4)                                                                                                                                                                                                                                
    for axis in ['top','bottom','left','right']:                                                                                                                                                                                                                  
      axinsets[-1].spines[axis].set_linewidth(0.3)                                                                                                                                                                                                                
    axinsets[-1].set_xticks([])
    boxoff(axinsets[-1])
    if (iax == 2 or iax == 3) and (iay == 2 or iay == 3):
      polygon = Polygon(array([[-2,4,4,-2],[0,0,60,60]]).T, True)
      p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
      p.set_facecolor('#EEEEEE')
      p.set_edgecolor('none')
      axinsets[-1].add_collection(p)
    if (iax == 4 or iax == 5) and (iay == 1):
      polygon = Polygon(array([[-2,4,4,-2],[0,0,60,60]]).T, True)
      p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
      p.set_facecolor('#FFE8E8')
      p.set_edgecolor('none')
      axinsets[-1].add_collection(p)

axarr = axs.reshape(prod(axs.shape),).tolist()


GluRCoeffs = ['0.5,0.5,1.5,1.5']
myfreqs = [0.001]+[1.0*i for i in range(1,21)]

min10 = [inf,inf,inf]
max10 = [-inf,-inf,-inf]
min20 = [inf,inf,inf]
max20 = [-inf,-inf,-inf]
contr10 = [-inf,-inf,-inf]
contr20 = [-inf,-inf,-inf]

blockeds = ['Calbin,CalbinC,DAGK,Gi,MGluR,NCX,Ng,PDE4,PKC,PP1', 'Calbin,CalbinC,DAGK,Gi,MGluR,NCX,Ng,PDE4,PKC,PP1', 'CK,DAGK,Gqabg,PKA', 'DAGK,PKA,PKC,PLA2,PP1', 'CK,Calbin,CalbinC,DAGK,Gi,Gqabg,MGluR,NCX,Ng,PDE4,PKA,PKC,PP1', 'Calbin,CalbinC,DAGK,Gi,MGluR,NCX,Ng,PDE4,PKA,PKC,PLA2,PP1', 'Calbin,CalbinC,DAGK,Gi,Ng,PDE4,PKC', 'DAGK,Gi,Ng,PDE4,PP1', 'CK,DAGK,Gqabg,PKA', 'DAGK,PKA,PKC,PLA2,PP1', 'CK,Calbin,CalbinC,DAGK,Gi,Gqabg,Ng,PDE4,PKA,PKC', 'DAGK,Gi,Ng,PDE4,PKA,PKC,PLA2,PP1', 'DAGK,MGluR,NCX,PDE4,PP1', 'MGluR', 'CK,DAGK,Gqabg,PKA', 'DAGK,PKA,PKC,PLA2,PP1', 'CK,DAGK,Gqabg,MGluR,NCX,PDE4,PKA,PP1', 'DAGK,MGluR,PKA,PKC,PLA2,PP1', 'DAGK,PDE4', 'CK,DAGK,Gqabg,PKA', 'DAGK,PKA,PKC,PLA2,PP1', 'CK,DAGK,Gqabg,PDE4,PKA', 'DAGK,PKA,PKC,PLA2,PP1']
blockedcoeffs = [[0.967,0.967,0.974,0.979,1.062,1.090,0.973,1.054,0.994,1.066], [1.025,1.025,0.988,0.968,1.093,1.042,0.941,0.967,0.984,0.986], [0.905,1.140,1.109,1.079], [1.158,1.055,1.056,0.864,0.897], [0.905,0.967,0.967,1.029,0.979,1.109,1.062,1.090,0.973,1.054,1.079,0.994,1.066], [1.025,1.025,1.045,0.968,1.093,1.042,0.941,0.967,1.055,1.020,0.864,0.942], [0.985,0.985,0.924,0.971,0.962,1.012,0.996], [0.992,0.984,0.944,0.977,0.992], [0.959,1.072,1.063,1.048], [1.134,1.039,1.037,0.933,0.940], [0.959,0.985,0.985,0.998,0.971,1.063,0.962,1.012,1.048,0.996], [1.039,0.984,0.944,0.977,1.039,1.037,0.933,0.966], [0.918,1.062,1.090,1.054,1.066], [1.093], [0.905,1.140,1.109,1.079], [1.158,1.055,1.056,0.864,0.897], [0.905,1.029,1.109,1.062,1.090,1.054,1.079,1.066], [1.158,1.093,1.055,1.056,0.864,0.897], [0.924,1.012], [0.959,1.072,1.063,1.048], [1.134,1.039,1.037,0.933,0.940], [0.959,0.998,1.063,1.012,1.048], [1.134,1.039,1.037,0.933,0.940]]
geneparams = ['sig0,imp0\n,GWAS,PFC','sig0,imp0\n,GWAS,ACC','sig0,imp0\n,CommonMind,PFC','sig0,imp0\n,CommonMind,ACC','sig0,imp0\n,Both,PFC','sig0,imp0\n,Both,ACC','sig0,imp1,\nGWAS,PFC','sig0,imp1,\nGWAS,ACC','sig0,imp1,\nCommonMind,PFC','sig0,imp1,\nCommonMind,ACC','sig0,imp1,\nBoth,PFC','sig0,imp1,\nBoth,ACC','sig1,imp0\n,GWAS,PFC','sig1,imp0\n,GWAS,ACC','sig1,imp0\n,CommonMind,PFC','sig1,imp0\n,CommonMind,ACC','sig1,imp0\n,Both,PFC','sig1,imp0\n,Both,ACC','sig1,imp1,\nGWAS,PFC','sig1,imp1,\nCommonMind,PFC','sig1,imp1,\nCommonMind,ACC','sig1,imp1,\nBoth,PFC','sig1,imp1,\nBoth,ACC']
blockedTitles = ['comb'+str(i) for i in range(0,23)]

Caflux = 150.0
for iglur in range(0,len(GluRCoeffs)):
  DATAS_0_CONTROL = ltdltpcurves_many(['GluR1,GluR1_memb,GluR2,GluR2_memb,Cax'+GluRCoeffs[iglur]+',1.0'], T, Caflux, Lflux, Gluflux, [''], freqs = myfreqs, dodraw = False)
  #print(str(DATAS_0_CONTROL))
  minmaxes = []
  
  for iblocked in range(0,len(blockeds)):
    iblockedCoeff = 1 if 'PFC' in geneparams[iblocked] else 2
    for iicoeff in range(0,1):
    #blocked_base = 'GluR1,GluR1_memb,GluR2,GluR2_memb,PDE1A,Cax0.5,0.5,1.5,1.5,0.8,1.0'
      #iblockedCoeff = blockedSet[iicoeff]
      #coeff = blockedCoeffs[iblockedCoeff]
      #Nsamespecies = blockeds[iblocked].count(',')+1
      blocked_base = 'GluR1,GluR1_memb,GluR2,GluR2_memb,CMcombx'+GluRCoeffs[iglur]+','+str(iblocked)
      DATAS_0 = ltdltpcurves_many([blocked_base], T, Caflux, Lflux, Gluflux, [''], freqs = myfreqs, dodraw = False)
      #print(str(DATAS_0))
      if DATAS_0[0][10][1] < min10[iglur]:
       min10[iglur] = DATAS_0[0][10][1]
      if DATAS_0[0][10][1] > max10[iglur]:
       max10[iglur] = DATAS_0[0][10][1]
      if DATAS_0[0][20][1] < min20[iglur]:
       min20[iglur] = DATAS_0[0][20][1]
      if DATAS_0[0][20][1] > max20[iglur]:
       max20[iglur] = DATAS_0[0][20][1]
      axarr[iblocked+(iblocked > 18)].plot(myfreqs,[DATAS_0[0][i][1] for i in range(0,len(DATAS_0[0]))],'b.-', lw=0.3, ms=1.0, mew=1.0, color=cols[iblockedCoeff], label = blockeds[iblocked])
      minmaxes.append([min([DATAS_0[0][i][1] for i in range(0,len(DATAS_0[0]))]),max([DATAS_0[0][i][1] for i in range(0,len(DATAS_0[0]))])])
      #if iglur == 0:
      axinsets[iblocked+(iblocked > 18)].bar(-1+2*iicoeff if iblocked < len(blockeds)-1 else iicoeff,DATAS_0[1][0][-1],color=cols[iblockedCoeff])
      #  #print("DATAS_0[0][0][0] = "+str(DATAS_0[0][0][0])+", DATAS_0[0][0][1] = "+str(DATAS_0[0][0][1]))
    axarr[iblocked+(iblocked > 18)].plot(myfreqs,[DATAS_0_CONTROL[0][i][1] for i in range(0,len(DATAS_0_CONTROL[0]))],'b.-', lw=0.3, ms=1.0, mew=1.0, color='#000000', label = 'Control')
    axarr[iblocked+(iblocked > 18)].set_xlim([0,20])
    axarr[iblocked+(iblocked > 18)].set_ylim([0.6,2.6])
    axinsets[iblocked+(iblocked > 18)].bar(0 if iblocked < len(blockeds)-1 else -1,DATAS_0_CONTROL[1][0][-1],color='#000000')
    axinsets[iblocked+(iblocked > 18)].set_xlim([-1.5,1.5])
    axinsets[iblocked+(iblocked > 18)].set_ylim([0,43])

    contr10[iglur] = DATAS_0_CONTROL[0][10][1]
    contr20[iglur] = DATAS_0_CONTROL[0][20][1]
    #axarr[iblocked].set_title(blockedTitles[iblocked],fontsize=4)
    mytxt = geneparams[iblocked][-3:]+', '+('Imputed' if 'imp1' in geneparams[iblocked] else 'Non-imputed')+',\n'+(geneparams[iblocked].replace('\n','').split(',')[2] if not 'Both' in geneparams[iblocked] else 'GWAS+CommonMind')+',\n'+\
            ('Only significant' if 'sig1' in geneparams[iblocked] else 'Significant +\nnon-significant')
    axarr[iblocked+(iblocked > 18)].text(0.5,max([DATAS_0_CONTROL[0][i][1] for i in range(0,len(DATAS_0_CONTROL[0]))])*1.02,mytxt,fontsize=4,ha='left',va='top')

  #axarr[1,iglur].set_xlabel(GluRCoeffs[iglur],fontsize=7)
print(str(min10)+","+str(max10)+","+str(min20)+","+str(max20))
print("control: "+str(contr10)+",2"+str(contr20))
print("minmaxes = "+str(minmaxes))
#axarr[0,1].legend(fontsize=5)
for iax in [0,1,2,3,4,6,7,8,9,10,12,13,14,15,16,20,21,22]:
  axarr[iax].set_xticks([0,10])
axarr[19].set_visible(False)
axinsets[19].set_visible(False)

f.savefig('fig_ltdltpcurves_CMcombs.eps')

