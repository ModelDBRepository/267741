#cp collectlocalcurrs_manyinputs_varyNMDA_oneISI_varyfreq.py collectlocalcurrsmut.py
#copied from  drawcurrs_stim_spines_info.py (which was copied from drawcurrs_stim_spines.py) 6.3.2019
import matplotlib
matplotlib.use('Agg')
from pylab import *
import scipy.io
from os.path import exists
#noisy_icell5_n100_5.0_Nsyn12_Econ0.0003_rateE0.7_seed17.mat
import mytools
import time

icell=0
mutID=0
stimN=100
stimfreq=1.0
Ninputs=1
Nsyn=50
Econ=0.0003
wNMDA=1.0
gNap=0.0
pulseamp = 1.0
dendtree = 'apic'
spinelocations = '100-200'
Npulses=1
neckLen=0.5
neckDiam=0.1

if len(sys.argv) > 1:
  icell = int(float(sys.argv[1]))
if len(sys.argv) > 2:
  mutID=int(float(sys.argv[2]))
if len(sys.argv) > 3:
  Nsyn = int(float(sys.argv[3]))
if len(sys.argv) > 4:
  Econ = float(sys.argv[4])
if len(sys.argv) > 5:
  wNMDA = float(sys.argv[5])
if len(sys.argv) > 6:
  gNap = float(sys.argv[6])
if len(sys.argv) > 7:
  pulseamp = float(sys.argv[7])
if len(sys.argv) > 8:
  dendtree = sys.argv[8]
if len(sys.argv) > 9:
  spinelocations = sys.argv[9]
if len(sys.argv) > 10:
  stimfreq = float(sys.argv[10])
if len(sys.argv) > 12:
  Npulses = int(float(sys.argv[12]))
if len(sys.argv) > 13:
  neckLen = float(sys.argv[13])
if len(sys.argv) > 14:
  neckDiam = float(sys.argv[14])
rateE=0.7
dtpulses=10.0
#pulseamps = [0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
#noisy_icell5_n100_5.0_neckLen0.5_neckDiam0.1_Nsyn5_Econ0.0003_rateE0.7_Npulses4_ISI60.0_dtpulses10.0_pulseamp1.0_seed13.mat
#ISIS=(-90.0 -80.0 -70.0 -60.0 -50.0 -40.0 -30.0 -20.0 -10.0 0.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0)
#ISIs = [-180.0, -160.0, -140.0, -120.0, -100.0, -80.0, -60.0, -50.0, -40.0, -30.0, -25.0, -20.0, -15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0]
ISIs = [-200.0, -180.0, -160.0, -140.0, -120.0, -100.0, -80.0, -60.0, -50.0, -40.0, -30.0, -25.0, -20.0, -15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0,-210.0, -220.0, -230.0, -2.5, 2.5]
cols = mytools.colorsredtolila(len(ISIs)+1,0.7)
t_index_start = 900 # tsamp_start is 4500, stimulus starts at 4500. By setting this 900, we will save clips from -100 to 900
T = int(1000/stimfreq)
if True:
    for iISI in range(0,len(ISIs)):
      if len(sys.argv) > 11 and float(sys.argv[11]) != ISIs[iISI]:
        continue
      dists = []
      currs = []
      currClips = []
      ISIpairing = ISIs[iISI]
      foundone = False
      for rdSeed in range(0,401):
        # #('noisy_icell'+str(icell)+'_n100_5.0_Nsyn12_Econ0.0003_rateE0.7_seed'+str(rdSeed)+'.mat'):
        if exists("noisy_icell"+str(icell)+"_imutc"+str(mutID)+"_n"+str(stimN)+"_"+str(stimfreq)+"_neckLen"+str(neckLen)+"_neckDiam"+str(neckDiam)+"_Nsyn"+str(Nsyn)+"_Ninputs"+str(Ninputs)+dendtree+spinelocations+"_Econ"+str(Econ)+'_wNMDA'+str(wNMDA)+"_gNap"+str(gNap)+"_rateE"+str(rateE)+"_Npulses"+str(Npulses)+"_ISI"+str(ISIpairing)+"_dtpulses"+str(dtpulses)+"_pulseamp"+str(pulseamp)+"_seed"+str(rdSeed)+".mat"):
          print("loading noisy_icell"+str(icell)+"_imutc"+str(mutID)+"_n"+str(stimN)+"_"+str(stimfreq)+"_neckLen"+str(neckLen)+"_neckDiam"+str(neckDiam)+"_Nsyn"+str(Nsyn)+"_Ninputs"+str(Ninputs)+dendtree+spinelocations+"_Econ"+str(Econ)+'_wNMDA'+str(wNMDA)+"_gNap"+str(gNap)+"_rateE"+str(rateE)+"_Npulses"+str(Npulses)+"_ISI"+str(ISIpairing)+"_dtpulses"+str(dtpulses)+"_pulseamp"+str(pulseamp)+"_seed"+str(rdSeed)+".mat")
          foundone = True
          A = scipy.io.loadmat("noisy_icell"+str(icell)+"_imutc"+str(mutID)+"_n"+str(stimN)+"_"+str(stimfreq)+"_neckLen"+str(neckLen)+"_neckDiam"+str(neckDiam)+"_Nsyn"+str(Nsyn)+"_Ninputs"+str(Ninputs)+dendtree+spinelocations+"_Econ"+str(Econ)+'_wNMDA'+str(wNMDA)+"_gNap"+str(gNap)+"_rateE"+str(rateE)+"_Npulses"+str(Npulses)+"_ISI"+str(ISIpairing)+"_dtpulses"+str(dtpulses)+"_pulseamp"+str(pulseamp)+"_seed"+str(rdSeed)+".mat")
          for isyn in range(0,Nsyn):
            currThis = A['DATA'][0,:]*0.0
            for iinput in range(0,Ninputs):
              currThis = currThis + A['DATA'][isyn*Ninputs+iinput,:]
            for ipulse in range(0,int((A['Nsamp']-t_index_start)/T)):
              currClips.append(currThis[t_index_start+T*ipulse:t_index_start+T*(ipulse+1)])
            currs.append(currThis[:])
            if A['stimdistances'].shape[0] == 1 and len(A['stimdistances'].shape) > 1 and A['stimdistances'].shape[1] > 1:
              dists.append(A['stimdistances'][0][isyn])
            else:
              dists.append(A['stimdistances'][isyn])
      if foundone:
        print("found at least one seed. saved ts from "+str(-500+t_index_start)+" to "+str(-500+t_index_start+T))
        #print array(currClips[:])
        #print array(currClips[:]).shape
        #print len(currClips)
        #print currClips[0].shape
        #print A['DATA'][0,:].shape
        print(A['DATA'][isyn,t_index_start+T*(ipulse-1):t_index_start+T*(ipulse-1+1)].shape)
        print(A['DATA'][isyn,t_index_start+T*ipulse:t_index_start+T*(ipulse+1)].shape)
        scipy.io.savemat('currClips'+str(icell)+"_imutc"+str(mutID)+'_neckLen'+str(neckLen)+'_neckDiam'+str(neckDiam)+'_stimfreq'+str(stimfreq)+'_pulseamp'+str(pulseamp)+'_Nsyn'+str(Nsyn)+'_Ninputs'+str(Ninputs)+dendtree+spinelocations+'_Econ'+str(Econ)+'_wNMDA'+str(wNMDA)+"_gNap"+str(gNap)+'_Npulses'+str(Npulses)+'_ISI'+str(ISIpairing)+'.mat',
                         {'currClips': mean(array(currClips[:]),axis=0), 'ts': list(range(-500+t_index_start,-500+t_index_start+T)), 'ISI': ISIpairing, 'filenameexample': "noisy_icell"+str(icell)+"_imutc"+str(mutID)+"_n"+str(stimN)+"_"+str(stimfreq)+"_neckLen"+str(neckLen)+"_neckDiam"+str(neckDiam)+"_Nsyn"+str(Nsyn)+"_Ninputs"+str(Ninputs)+dendtree+spinelocations+"_Econ"+str(Econ)+'_wNMDA'+str(wNMDA)+"_gNap"+str(gNap)+"_rateE"+str(rateE)+"_Npulses"+str(Npulses)+"_ISI"+str(ISIpairing)+"_dtpulses"+str(dtpulses)+"_pulseamp"+str(pulseamp)+"_seed"+str(rdSeed)+".mat"})
        print('to be removed: rm '+"noisy_icell"+str(icell)+"_imutc"+str(mutID)+"_n"+str(stimN)+"_"+str(stimfreq)+"_neckLen"+str(neckLen)+"_neckDiam"+str(neckDiam)+"_Nsyn"+str(Nsyn)+"_Ninputs"+str(Ninputs)+dendtree+spinelocations+"_Econ"+str(Econ)+'_wNMDA'+str(wNMDA)+"_gNap"+str(gNap)+"_rateE"+str(rateE)+"_Npulses"+str(Npulses)+"_ISI"+str(ISIpairing)+"_dtpulses"+str(dtpulses)+"_pulseamp"+str(pulseamp)+"_seed*.mat")
      else:
        print("noisy_icell"+str(icell)+"_imutc"+str(mutID)+"_n"+str(stimN)+"_"+str(stimfreq)+"_neckLen"+str(neckLen)+"_neckDiam"+str(neckDiam)+"_Nsyn"+str(Nsyn)+"_Ninputs"+str(Ninputs)+dendtree+spinelocations+"_Econ"+str(Econ)+'_wNMDA'+str(wNMDA)+"_gNap"+str(gNap)+"_rateE"+str(rateE)+"_Npulses"+str(Npulses)+"_ISI"+str(ISIpairing)+"_dtpulses"+str(dtpulses)+"_pulseamp"+str(pulseamp)+"_seed*.mat not found")
      #currClips_all.append(array(currClips[:]))
      #currClipmeans_all.append(mean(array(currClips[:]),axis=0))
      

