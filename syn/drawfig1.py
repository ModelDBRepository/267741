from pylab import *
import scipy.io
from os.path import exists
import os
import mytools
import calcconds

#ltdltpcurve: get the conductances, GluR1 concentrations at the membrane and active PKA and PKC concentrations for all stimulation frequencies
def ltdltpcurve(blocked='GluR1,GluR1_memb,GluR2,GluR2_memb,Cax0.5,0.5,1.5,1.5,1.0', T = 100, Caflux = 150.0, Lflux = 5.0, Gluflux = 10.0, altered = 'k1x1.0', freqs = [0.5, 1.0, 2.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 22.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 200.0, 300.0, 500.0],dodraw = True):

 Nstims = [int(T*x) for x in freqs]

 tposts = [20, 30, 40] #min after stimulus onset

 condposts = []
 condposts_abs = []
 GluR1posts = []
 GluR1S831posts = []
 GluR2posts = []
 PKCposts = []
 PKAposts = []
 GluR1homposts = []

 addition = ''
 for ifreq in range(0,len(freqs)):
  freq = freqs[ifreq]
  Nstim = Nstims[ifreq]
  filename = addition+'nrn_tstop27000000_tol1e-06_'+blocked+'_'+altered+'_onset24040000.0_n'+str(Nstim)+'_freq'+str(freq)+'_dur3.0_flux'+str(Caflux)+'_Lflux'+str(Lflux)+'_Gluflux'+str(Gluflux)+'_AChflux'+str(Gluflux)+'_Ntrains1_trainT1000.0.mat'
  if ifreq == 0 and not exists(filename) and not exists(filename+'.mat') and (exists('mats/'+filename) or exists('mats/'+filename+'.mat')): #do this only once to save i/o, the same should apply to either all freqs or none
    addition = 'mats/'
    filename = addition+'nrn_tstop27000000_tol1e-06_'+blocked+'_'+altered+'_onset24040000.0_n'+str(Nstim)+'_freq'+str(freq)+'_dur3.0_flux'+str(Caflux)+'_Lflux'+str(Lflux)+'_Gluflux'+str(Gluflux)+'_AChflux'+str(Gluflux)+'_Ntrains1_trainT1000.0.mat'
  if 'control' in altered:
    filename_alt = 'nrn_'+str(Nstim)+'_'+str(freq)+'_'+str(Caflux)+'_'+str(Lflux)+'_'+str(Gluflux)+'_'+str(Gluflux)+'_'+blocked+'_'+altered+'.mat'
  else:  
    filename_alt = 'nrn_altered2_'+str(Nstim)+'_'+str(freq)+'_'+str(Caflux)+'_'+str(Lflux)+'_'+str(Gluflux)+'_'+str(Gluflux)+'_'+blocked+'_'+altered+'.mat'
  if not exists(filename) and exists(filename+'.mat'):
    filename = filename+'.mat'
  elif not exists(filename):
    filename = filename_alt
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
    iPKA = [i for i in range(0,len(A['headers'])) if 'PKAc' in A['headers'][i] and 'cAMP' not in A['headers'][i]]
    iPKC = [i for i in range(0,len(A['headers'])) if 'PKCt' in A['headers'][i]]
    GluR1s = [sum([A['DATA'][i][it] for i in iGluR1]) for it in range(0,len(times))]
    GluR1S831s = [sum([A['DATA'][i][it] for i in iGluR1S831]) for it in range(0,len(times))]
    GluR2s = [sum([A['DATA'][i][it] for i in iGluR2]) for it in range(0,len(times))]
    PKAs = [sum([A['DATA'][i][it] for i in iPKA]) for it in range(0,len(times))]
    PKCs = [sum([A['DATA'][i][it] for i in iPKC]) for it in range(0,len(times))]
    GluR1homs = [(GluR1s[i]+GluR2s[i])/4.0 * GluR1s[i]**4/(GluR1s[i]+GluR2s[i])**4 for i in range(0,len(times))]

    GluR1posts.append([GluR1s[its[ipost]] for ipost in range(0,len(tposts))])
    GluR1S831posts.append([GluR1S831s[its[ipost]] for ipost in range(0,len(tposts))])
    GluR2posts.append([GluR2s[its[ipost]] for ipost in range(0,len(tposts))])
    PKAposts.append([max(PKAs) for ipost in range(0,len(tposts))])
    PKCposts.append([max(PKCs) for ipost in range(0,len(tposts))])
    GluR1homposts.append([GluR1homs[its[ipost]] for ipost in range(0,len(tposts))])
  else:
    print(filename+' not found')
    condposts.append(list(nan*ones([len(tposts),])))
    condposts_abs.append(list(nan*ones([len(tposts),])))
    GluR1posts.append(list(nan*ones([len(tposts),])))
    GluR1S831posts.append(list(nan*ones([len(tposts),])))
    GluR2posts.append(list(nan*ones([len(tposts),])))
    PKAposts.append(list(nan*ones([len(tposts),])))
    PKCposts.append(list(nan*ones([len(tposts),])))
    GluR1homposts.append(list(nan*ones([len(tposts),])))

 return [condposts, condposts_abs, GluR1posts, GluR2posts, GluR1S831posts, GluR1homposts, PKAposts, PKCposts]

#ltdltpcurves_many: get the LTD/LTP curve for many different versions of the model (altereds can have more than one case)
def ltdltpcurves_many(blockeds=['GluR1,GluR1_memb,GluR2,GluR2_memb,Cax0.5,0.5,1.5,1.5,1.0'], T = 100, Caflux = 150.0, Lflux = 5.0, Gluflux = 10.0, altereds = ['k1x1.0'], freqs = [0.5, 1.0, 2.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 22.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 200.0, 300.0, 500.0],dodraw = True, epsfilename = '', givenlabels = []):
 tposts = [20, 30, 40] #min after stimulus onset
 if len(epsfilename) == 0:
   epsfilename = 'ltdltpcurve_varyNstim_many_'+blockeds[0]+'_'+altereds[0]+'_'+str(T)+'_Caflux'+str(Caflux)+'_Lflux'+str(Lflux)+'_Gluflux'+str(Gluflux)+'.eps'
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



f,axarr = subplots(9,1)
axarr[0].set_position([0.06,0.49,0.32,0.46])
for ix in range(0,3):
  for iy in range(0,2):
    axarr[1+3*iy+ix].set_position([0.44+0.20*ix,0.49+0.26*(1-iy),0.14,0.2])
axarr[7].set_position([0.12+0.29*0,0.41+0.26*(-1),0.14,0.2])
axarr[8].set_position([0.12+0.29*1,0.41+0.26*(-1),0.14,0.2])
axnew = []
axnew.append(f.add_axes([0.118, 0.77, 0.085, 0.1]))
axnew.append(f.add_axes([0.288, 0.57, 0.085, 0.1]))
axnew.append(f.add_axes([0.44+0.20*2+0.08, 0.77+0.1, 0.05, 0.06]))
axnew.append(f.add_axes([0.44+0.20*2+0.08, 0.51+0.07, 0.05, 0.06]))

for ax in list(axarr) + axnew:
  for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(4)
    for axis in ['top','bottom','left','right']:
      ax.spines[axis].set_linewidth(0.5)

for ax in axnew:
  for q in ['top','right']:
    ax.spines[q].set_visible(False)

T = 100
    
altereds_all = ['control','icase0','icase1','icase2','icase3']
myfreqs = [0.001]+[1.0*i for i in range(1,28)]

casenames = ['control','no PKA bind GluR1','no PKA bind I1','no PKA bind GluR1 or I1','no PKC bind GluR2','no PKC bind GluR1 or GluR2','no CK bind GluR1','no PKC bind GluR1','no CK or PKC bind GluR1','control']

Caflux = 150.0
Lflux = 5.0
Gluflux = 10.0
blocked_base = 'GluR1,GluR1_memb,GluR2,GluR2_memb,Ca,Cax0.5,0.5,1.5,1.5,1.0,1.0'
DATAS_0 = ltdltpcurves_many([blocked_base], T, Caflux, Lflux, Gluflux, altereds_all, freqs = myfreqs, dodraw = False, epsfilename = 'ltdltpcurves2_altereds_'+str(T)+'_'+str(Caflux)+'_'+str(Lflux)+'_'+str(Gluflux)+'_'+blocked_base+'.eps',givenlabels = casenames)
axarr[0].plot(myfreqs,[DATAS_0[0][0][i][1] for i in range(0,len(DATAS_0[0][0]))],'b.-')
axarr[0].plot([myfreqs[i] for i in [5,15]],[DATAS_0[0][0][i][1] for i in [5,15]],'b.',ms=9)
axarr[0].set_title('Relative conductance',fontsize=7)
axarr[0].set_ylabel('fold change',fontsize=7)

axarr[1].plot(myfreqs,[DATAS_0[4][0][i][1] for i in range(0,len(DATAS_0[0][0]))],'cx-',ms=1.5,mew=1.5,lw=0.3)
axarr[1].plot(myfreqs,[DATAS_0[0][0][i][1] for i in range(0,len(DATAS_0[0][0]))],'b.-',ms=1.5,mew=1.5,lw=0.3)
axarr[1].set_title('Rel. conductance',fontsize=7)
axarr[1].set_ylabel('fold change',fontsize=7)

axarr[2].plot(myfreqs,[1000*DATAS_0[4][7][i][1] for i in range(0,len(DATAS_0[0][0]))],'cx-',ms=1.5,mew=1.5,lw=0.3)
axarr[2].plot(myfreqs,[1000*DATAS_0[0][7][i][1] for i in range(0,len(DATAS_0[0][0]))],'b.-',ms=1.5,mew=1.5,lw=0.3)
axarr[2].set_title('Active PKC',fontsize=7)
axarr[2].set_ylabel('$\mu$M',fontsize=7)

axarr[3].plot(myfreqs,[1e6*DATAS_0[4][3][i][1] for i in range(0,len(DATAS_0[0][0]))],'cx-',ms=1.5,mew=1.5,lw=0.3)
axarr[3].plot(myfreqs,[1e6*DATAS_0[0][3][i][1] for i in range(0,len(DATAS_0[0][0]))],'b.-',ms=1.5,mew=1.5,lw=0.3)
axarr[3].set_title('GluR2 at memb.',fontsize=7)
axarr[3].set_ylabel('nM',fontsize=7)

axarr[4].plot(myfreqs,[DATAS_0[1][0][i][1] for i in range(0,len(DATAS_0[0][0]))],'rx-',ms=1.5,mew=1.5,lw=0.3)
axarr[4].plot(myfreqs,[DATAS_0[0][0][i][1] for i in range(0,len(DATAS_0[0][0]))],'b.-',ms=1.5,mew=1.5,lw=0.3)
axarr[4].set_title('Rel. conductance',fontsize=7)
axarr[4].set_ylabel('fold change',fontsize=7)

axarr[5].plot(myfreqs,[1000*DATAS_0[1][6][i][1] for i in range(0,len(DATAS_0[0][0]))],'rx-',ms=1.5,mew=1.5,lw=0.3)
axarr[5].plot(myfreqs,[1000*DATAS_0[0][6][i][1] for i in range(0,len(DATAS_0[0][0]))],'b.-',ms=1.5,mew=1.5,lw=0.3)
axarr[5].set_title('Active PKA',fontsize=7)
axarr[5].set_ylabel('$\mu$M',fontsize=7)

axarr[6].plot(myfreqs,[1e6*DATAS_0[1][2][i][1] for i in range(0,len(DATAS_0[0][0]))],'rx-',ms=1.5,mew=1.5,lw=0.3)
axarr[6].plot(myfreqs,[1e6*DATAS_0[0][2][i][1] for i in range(0,len(DATAS_0[0][0]))],'b.-',ms=1.5,mew=1.5,lw=0.3)
axarr[6].set_title('GluR1 at memb.',fontsize=7)
axarr[6].set_ylabel('nM',fontsize=7)

Lflux = 0.2
blocked_base = 'GluR1,GluR1_memb,GluR2,GluR2_memb,PP1,Cax2.0,2.0,0.0,0.0,0.5,1.0'
DATAS_1 = ltdltpcurves_many([blocked_base], T, Caflux, Lflux, Gluflux, altereds_all, freqs = myfreqs, dodraw = False, epsfilename = 'ltdltpcurves2_altereds_'+str(T)+'_'+str(Caflux)+'_'+str(Lflux)+'_'+str(Gluflux)+'_'+blocked_base+'.eps',givenlabels = casenames)
axarr[7].plot(myfreqs,[DATAS_1[0][0][i][1] for i in range(0,len(DATAS_1[0][0]))],'b.-',ms=1.5,mew=1.5,lw=0.3)
axarr[7].set_title('Rel. conductance (R1 only, 50% less\nPP1, milder beta-adr. stimulation)',fontsize=7)
axarr[7].set_ylabel('fold change',fontsize=7)

Lflux = 5.0
blocked_base = 'GluR1,GluR1_memb,GluR2,GluR2_memb,R,PP1,Cax2.0,2.0,0.0,0.0,0.0,0.0,1.0'
DATAS_2 = ltdltpcurves_many([blocked_base], T, Caflux, Lflux, Gluflux, altereds_all, freqs = myfreqs, dodraw = False, epsfilename = 'ltdltpcurves2_altereds_'+str(T)+'_'+str(Caflux)+'_'+str(Lflux)+'_'+str(Gluflux)+'_'+blocked_base+'.eps',givenlabels = casenames)
axarr[8].plot(myfreqs,[DATAS_2[0][0][i][1] for i in range(0,len(DATAS_2[0][0]))],'b.-',ms=1.5,mew=1.5,lw=0.3)
axarr[8].set_title('Rel. conductance (R1 only,\nno $\\beta$-AR activation, no PP1)',fontsize=7)
axarr[8].set_ylabel('fold change',fontsize=7)

#f.text(0.6, 0.35, 'control', fontsize=8)
myleg = mytools.mylegend(f,[0.6,0.25,0.3,0.12],['b-','c-','r-'],['control','no S880 phosporylation by PKC','no S845 phosphorylation by PKA'],1,2,0.5,0.35,myfontsize=8)
for q in ['top','right','bottom','left']:
  myleg.spines[q].set_visible(False)

for iax in range(0,7):
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.03, pos.y1, chr(ord('A')+iax), fontsize=12)
for iax in range(7,9):
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.1+0.04*(iax==8), pos.y1+0.02, chr(ord('A')+iax), fontsize=12)

for iax in [0]+list(range(4,9)):
  axarr[iax].set_xlabel('$f$ (Hz)', fontsize=7)

conds,times = calcconds.calcconds_nrn('nrn_500_5.0_100.0_5.0_10.0_10.0_GluR1,GluR1_memb,GluR2,GluR2_memb,Ca,Cax0.5,0.5,1.5,1.5,1.0,1.0_control.mat')
axnew[0].plot([(x-24040000)/60000 for x in times],conds,'b-',lw=0.5)
conds,times = calcconds.calcconds_nrn('nrn_1500_15.0_100.0_5.0_10.0_10.0_GluR1,GluR1_memb,GluR2,GluR2_memb,Ca,Cax0.5,0.5,1.5,1.5,1.0,1.0_control.mat')
axnew[1].plot([(x-24040000)/60000 for x in times],conds,'b-',lw=0.5)

for iax in range(0,2):
  axnew[iax].set_xlim([0,40])
  axnew[iax].set_xticks([0,15,30])
  axnew[iax].set_xlabel('$t$ (min)',fontsize=7)
  axnew[iax].set_ylabel('(pS)',fontsize=7)
axnew[0].text(15,34,'$f$ = 5 Hz',fontsize=5,ha='left',va='center')
axnew[1].text(15,34,'$f$ = 15 Hz',fontsize=5,ha='left',va='center')

axnew[2].plot(myfreqs, [(1e6*DATAS_0[0][2][i][1])**4/(1e6*DATAS_0[0][2][i][1]+1e6*DATAS_0[0][3][i][1])**4 for i in range(0,len(DATAS_0[0][0]))],'b.-',ms=1.0,mew=1.0,lw=0.3)
axnew[2].plot(myfreqs, [(1e6*DATAS_0[4][2][i][1])**4/(1e6*DATAS_0[4][2][i][1]+1e6*DATAS_0[4][3][i][1])**4 for i in range(0,len(DATAS_0[0][0]))],'cx-',ms=1.0,mew=1.0,lw=0.3)

axnew[3].plot(myfreqs, [(1e6*DATAS_0[0][2][i][1])**4/(1e6*DATAS_0[0][2][i][1]+1e6*DATAS_0[0][3][i][1])**4 for i in range(0,len(DATAS_0[0][0]))],'b.-',ms=1.0,mew=1.0,lw=0.3)
axnew[3].plot(myfreqs, [(1e6*DATAS_0[1][2][i][1])**4/(1e6*DATAS_0[1][2][i][1]+1e6*DATAS_0[1][3][i][1])**4 for i in range(0,len(DATAS_0[0][0]))],'rx-',ms=1.0,mew=1.0,lw=0.3)

for iax in range(2,4):
  axnew[iax].set_xlim([0,20])
  axnew[iax].set_xticks([0,20])
  axnew[iax].set_xlabel('$f$ (Hz)',fontsize=7)
  axnew[iax].set_ylim([0,0.2])

for ax in axarr:
 ax.set_xlim([0,20])

axnew[3].set_visible(False) #This is not needed

f.savefig('fig_ltdltpcurves.eps')

