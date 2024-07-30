from pylab import *
import scipy.io
from sklearn.linear_model import LinearRegression
import pandas as pd
from os.path import exists

f,axarr = subplots(1,3)
for iax in range(0,3):
  if iax < 2:
    axarr[iax].set_position([0.08+0.24*iax,0.45,0.17,0.48])
    axarr[iax].set_ylim([-0.13,0.18])
    axarr[iax].plot([-0.5,5.5],[0,0],'k-',lw=0.5)
    axarr[iax].spines['top'].set_visible(False)
    axarr[iax].spines['right'].set_visible(False)
    axarr[iax].spines['bottom'].set_visible(False)
    axarr[iax].get_xaxis().tick_bottom()
    axarr[iax].get_yaxis().tick_left()
    axarr[iax].set_xticks([])
  else:
    axarr[iax].set_position([0.08+0.24*iax,0.45,0.4,0.48])
      
  for tick in axarr[iax].xaxis.get_major_ticks() + axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(4)

  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.05, pos.y1 - 0.01, chr(ord('A')+iax), fontsize=10)
  

A = scipy.io.loadmat('corrs_lips.mat')
thrs = [1e-5, 5e-6, 1e-6]


PRSlabels_all = ['genes_lips_synaptic_ABFGJLOQ']

PRSlabel_names = ['genes_new', 'genes_new_synaptic', 'genes_all_synaptic_PKA', 'genes_all_synaptic_PKC', 'genes_new_ionchannels'] 
iPRSlabels = [[i for i in range(0,len(PRSlabels_all)) if PRSlabels_all[i] == PRSlabel_names[j]][0] for j in range(0,len(PRSlabel_names))]
PRSlabels = ['Plasticity + ion channels', 'Plasticity', 'Plasticity, only PKA', 'Plasticity, only PKC', 'Ion channels']

corrs = A['corrs']
pvals = [[A['betas_pvals'][i,j][:,1,1] for j in range(0,A['betas_pvals'].shape[1])] for i in range(0,A['betas_pvals'].shape[0])]
betas = [[A['betas_pvals'][i,j][:,0,1] for j in range(0,A['betas_pvals'].shape[1])] for i in range(0,A['betas_pvals'].shape[0])]

cols = ['#6666AA','#999944','#229922','#BB3333','#993399','#999999','#777700','#999955']
ithrss = [0]
ielphys = 8
for iithrs in range(0,len(ithrss)):
  ithrs = ithrss[iithrs]
  for iiPRS in range(0,len(iPRSlabels)):
    print(str(pvals[ielphys][iPRSlabels[iiPRS]][ithrs]))
    if iithrs < 2:
      axarr[iithrs].bar(iiPRS,A['corrs'][ielphys,iPRSlabels[iiPRS]][0][ithrs],facecolor=cols[iiPRS])
      axarr[iithrs].text(iiPRS,max(0,A['corrs'][ielphys,iPRSlabels[iiPRS]][0][ithrs])+0.015,PRSlabels[iiPRS],fontsize=4.4,rotation=90,ha='center',va='bottom')
      fontweight = 'normal' if pvals[ielphys][iPRSlabels[iiPRS]][ithrs] >= 0.05 else 'bold'
      fontcolor = '#303030' if pvals[ielphys][iPRSlabels[iiPRS]][ithrs] >= 0.05 else '#000000'
      if abs(A['corrs'][ielphys,iPRSlabels[iiPRS]][0][ithrs]) > 0.045:
        if A['corrs'][ielphys,iPRSlabels[iiPRS]][0][ithrs] > 0:
          axarr[iithrs].text(iiPRS,0.004,'p='+'{:.3f}'.format(pvals[ielphys][iPRSlabels[iiPRS]][ithrs]),fontsize=4.4,rotation=90,ha='center',va='bottom',weight=fontweight,color=fontcolor)
        else:
          axarr[iithrs].text(iiPRS,-0.004,'p='+'{:.3f}'.format(pvals[ielphys][iPRSlabels[iiPRS]][ithrs]),fontsize=4.4,rotation=90,ha='center',va='top',weight=fontweight,color=fontcolor)
      else:
        if A['corrs'][ielphys,iPRSlabels[iiPRS]][0][ithrs] > 0:
          axarr[iithrs].text(iiPRS,-0.004,'p='+'{:.3f}'.format(pvals[ielphys][iPRSlabels[iiPRS]][ithrs]),fontsize=4.4,rotation=90,ha='center',va='top',weight=fontweight,color=fontcolor)
        else:
          axarr[iithrs].text(iiPRS,A['corrs'][ielphys,iPRSlabels[iiPRS]][0][ithrs]-0.004,'p='+'{:.3f}'.format(pvals[ielphys][iPRSlabels[iiPRS]][ithrs]),fontsize=4.4,rotation=90,ha='center',va='top',weight=fontweight,color=fontcolor)

PRSlabel_name = 'genes_lips_synaptic_ABFGJLOQ'
iPRSlabel = [i for i in range(0,len(PRSlabels_all)) if PRSlabels_all[i] == PRSlabel_name][0]
if exists('PRSs_'+PRSlabels_all[iPRSlabel]+'_elphys8_thr3.mat'):
  print('Loading PRSs_'+PRSlabels_all[iPRSlabel]+'_elphys8_thr3.mat')
  B = scipy.io.loadmat('PRSs_'+PRSlabels_all[iPRSlabel]+'_elphys8_thr3.mat')
  axarr[2].plot(B['PRS'][0],B['Y'][0],'k.',color=cols[5],lw=0.9,mew=0.9,ms=0.9)
  Xdf = pd.DataFrame(B['PRS'][0],columns = ['PRS'])
  Ydf = pd.DataFrame(B['Y'][0],columns = ['Y'])

  clf = LinearRegression().fit(Xdf, array(Ydf))

  print('Beta = '+str(clf.coef_)+" or "+str(betas[ielphys][iPRSlabel][ithrs])+", mycorr = "+str(A['corrs'][ielphys,iPRSlabel][0][ithrs])+", p = "+str(pvals[ielphys][iPRSlabel][ithrs])+", ielphys="+str(ielphys)+", "+PRSlabels_all[iPRSlabel])

  coefs = A['betas_pvals'][ielphys,iPRSlabel][ithrs,0,:]
  axarr[2].plot([-0.02,0.02],[coefs[0]-0.02*coefs[1], coefs[0]+0.02*coefs[1]],'k-',lw=0.3)

  axarr[2].set_title('Association between PRS(Synaptic, Lips et al.) and VEP-P1',fontsize=6)
  axarr[2].text(0.02,13,'corr='+'{:.2f}'.format(A['corrs'][ielphys,iPRSlabel][0][ithrs]),fontsize=6,ha='right')
  axarr[2].text(0.02,12,'p='+'{:.3f}'.format(pvals[ielphys][iPRSlabel][ithrs]),fontsize=6,ha='right',weight='bold')
else:
  print('PRSs_'+PRSlabels_all[iPRSlabel]+'_elphys8_thr3.mat does not exist')

axarr[0].set_ylabel('Correlation with VEP-P1',fontsize=6)
axarr[1].set_ylabel('Correlation with VEP-P1',fontsize=6)
axarr[2].set_ylabel('Correlation with VEP-P1',fontsize=6)

axarr[2].set_ylabel('VEP-P1 (A.U.)',fontsize=6)
axarr[2].set_xlabel('PRS',fontsize=6)

axarr[0].set_title('SNP inclusion\nthreshold $5\\times10^{-6}$',fontsize=6)
axarr[1].set_title('SNP inclusion\nthreshold $10^{-6}$',fontsize=6)
f.savefig('fig_correlations_PRS_EEGphenotype_suppl.eps')
