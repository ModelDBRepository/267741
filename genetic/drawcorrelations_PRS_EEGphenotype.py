#cp drawcorrelations_PRS_EEGphenotype.py drawcorrelations_PRS_EEGphenotype2.py
from pylab import *
import scipy.io
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import pandas as pd
from os.path import exists

f,axarr = subplots(1,9)
for iax in range(0,9):
  if iax < 3:
    axarr[iax].set_position([0.08+0.25*iax,0.48,0.15,0.472])
    axarr[iax].set_ylim([-0.13,0.16])
    axarr[iax].plot([-0.5,5.5],[0,0],'k-',lw=0.5)
    axarr[iax].spines['top'].set_visible(False)
    axarr[iax].spines['right'].set_visible(False)
    axarr[iax].spines['bottom'].set_visible(False)
    axarr[iax].get_xaxis().tick_bottom()
    axarr[iax].get_yaxis().tick_left()
    axarr[iax].set_xticks([])
  else:
    axarr[iax].set_position([0.07+0.12*(iax-3),0.26,0.07,0.16])
    axarr[iax].set_xlim([-0.02,0.02])
    axarr[iax].set_ylim([-14,14])
      
  for tick in axarr[iax].xaxis.get_major_ticks() + axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(4.8)

  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.035 - 0.015*(iax<3), pos.y1 - 0.007, chr(ord('A')+iax), fontsize=10)
  

A = scipy.io.loadmat('corrs.mat')
thrs = [1e-5, 5e-6, 1e-6]
ithrs = 0

PRSlabel_names = ['genes_new', 'genes_new_synaptic', 'genes_all_synaptic_PKA', 'genes_all_synaptic_PKC', 'genes_new_ionchannels'] 
PRSlabels = ['Plasticity + ion channels', 'Plasticity', 'Plasticity, only PKA', 'Plasticity, only PKC', 'Ion channels']

corrs = A['corrs']
pvals = [[A['betas_pvals'][i,j][:,1,1] for j in range(0,A['betas_pvals'].shape[1])] for i in range(0,A['betas_pvals'].shape[0])]
betas = [[A['betas_pvals'][i,j][:,0,1] for j in range(0,A['betas_pvals'].shape[1])] for i in range(0,A['betas_pvals'].shape[0])]

cols = ['#6666AA','#999944','#229922','#BB3333','#993399','#999999','#777700','#999955']
ielphys_plotted = [4,8,12]

for iielphys in range(0,len(ielphys_plotted)):
  ielphys = ielphys_plotted[iielphys]
  for iPRS in range(0,len(PRSlabels)):
    axarr[iielphys].bar(iPRS,A['corrs'][ielphys,iPRS][0][ithrs],facecolor=cols[iPRS])
    axarr[iielphys].text(iPRS,max(0,A['corrs'][ielphys,iPRS][0][ithrs])+0.015,PRSlabels[iPRS],fontsize=4.4,rotation=90,ha='center',va='bottom')
    print(str(pvals[ielphys][iPRS][ithrs]))
    fontweight = 'normal' if pvals[ielphys][iPRS][ithrs] >= 0.05 else 'bold'
    fontcolor = '#303030' if pvals[ielphys][iPRS][ithrs] >= 0.05 else '#000000'
    if abs(A['corrs'][ielphys,iPRS][0][ithrs]) > 0.045:
      if A['corrs'][ielphys,iPRS][0][ithrs] > 0:
        axarr[iielphys].text(iPRS,0.004,'p='+'{:.3f}'.format(pvals[ielphys][iPRS][ithrs]),fontsize=4.4,rotation=90,ha='center',va='bottom',weight=fontweight,color=fontcolor)
      else:
        axarr[iielphys].text(iPRS,-0.004,'p='+'{:.3f}'.format(pvals[ielphys][iPRS][ithrs]),fontsize=4.4,rotation=90,ha='center',va='top',weight=fontweight,color=fontcolor)
    else:
      if A['corrs'][ielphys,iPRS][0][ithrs] > 0:
        axarr[iielphys].text(iPRS,-0.004,'p='+'{:.3f}'.format(pvals[ielphys][iPRS][ithrs]),fontsize=4.4,rotation=90,ha='center',va='top',weight=fontweight,color=fontcolor)
      else:
        axarr[iielphys].text(iPRS,A['corrs'][ielphys,iPRS][0][ithrs]-0.004,'p='+'{:.3f}'.format(pvals[ielphys][iPRS][ithrs]),fontsize=4.4,rotation=90,ha='center',va='top',weight=fontweight,color=fontcolor)

    if ielphys == 8:
      if not exists('PRSs_'+PRSlabels_all[iPRS]+'_elphys8_thr3.mat'):
        print('PRSs_'+PRSlabels_all[iPRS]+'_elphys8_thr3.mat does not exist')
        continue
      B = scipy.io.loadmat('PRSs_'+PRSlabels_all[iPRS]+'_elphys8_thr3.mat')
      axarr[3+iPRS].plot(B['PRS'][0],B['Y'][0],'k.',color=cols[iPRS],lw=0.7,mew=0.7,ms=0.7)
      Xdf = pd.DataFrame(B['PRS'][0],columns = ['PRS'])
      Ydf = pd.DataFrame(B['Y'][0],columns = ['Y'])

      regr = LinearRegression()
      clf = regr.fit(Xdf, array(Ydf))
      r2 = r2_score(Ydf, regr.predict(Xdf))
      mycorr = corrcoef(array(Ydf)[:,0], regr.predict(Xdf)[:,0])[0,1]

      print('Beta = '+str(clf.coef_)+" or "+str(betas[ielphys][iPRS][ithrs])+", R2 = "+str(r2)+", mycorr = "+str(A['corrs'][ielphys,iPRS][0][ithrs])+", p = "+str(pvals[ielphys][iPRS][ithrs])+", mycorr(pred) = "+str(mycorr)+", iielphys="+str(iielphys)+", "+PRSlabels_all[iPRS])

      coefs = A['betas_pvals'][ielphys,iPRS][ithrs,0,:]
      axarr[3+iPRS].plot([-0.02,0.02],[coefs[0]-0.02*coefs[1], coefs[0]+0.02*coefs[1]],'k-',lw=0.3)

axarr[0].set_ylabel('Correlation with VEP-C1',fontsize=6)
axarr[1].set_ylabel('Correlation with VEP-P1',fontsize=6)
axarr[2].set_ylabel('Correlation with VEP-N1',fontsize=6)

axarr[3].set_ylabel('VEP-P1 (A.U.)',fontsize=6)
for iax in [3,4,5,6,7,8]:
  #axarr[iax].set_xlabel(('\n' if iax%2==0 else '')+'PRS('+PRSlabels[iax-3]+')',fontsize=6)
  axarr[iax].set_xlabel('PRS\n'+PRSlabels_linebreak[iax-3],fontsize=6)
f.savefig('fig_correlations_PRS_EEGphenotype2.pdf')
