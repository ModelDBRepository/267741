#changed from runmodelmut_old.py (11.8.2021): Use changes to total CaHVA, CaLVA, Ih, KCNB1 and SK conductances, as well as gamma_CaDyn
#cp runmodelsyns_exconly_stim_localspines_manyinputs_varyNMDA.py runmodelmut.py
#cp runmodelsyns_exconly_stim_spines_manyinputs_varyNMDA.py runmodelsyns_exconly_stim_localspines_manyinputs_varyNMDA.py
# - constrained where the spines can be
#difference between runmodelsyns_varEIsyn and runmodelsyns_varEIsynboth is that here their relation is kept fixed.
from neuron import h
import sys
from pylab import *
import time
import scipy.io
import mytools
import pickle

icell = 0 # Note: here icell=0 is the first cell, in ../l23pc2 icell=1 is the first cell
if len(sys.argv) > 1:
  icell = int(sys.argv[1])

tstop = 20000
rdSeed = 1
rateE = 0.0

mutID = -1
stimfreq = 1.0
stimN = 100
gCoeff = 1.0
Econ = 0.0003
weight_factor_NMDA = 1.0
gNap = 0.0
Nsyn = 5
Npulses = 4       #How many pulses in each burst paired (at perisomatic region) with presynaptic stimulation
dtpulses = 10.0   #Interburst interval in the perisomatic stimulus
pulsedur = 2.0    #How long is each current pulse
pulseamp = 1.0    #Amplitude in nanoamps
ISIpairing = 0.0  #The pairing (STDP) interval. 0 means that the first pulse starts at the same time as the presynaptic stimulus.

neckLen = 0.5 #https://books.google.no/books?id=reRPgWOImqIC&pg=PA149&lpg=PA149&dq=spine+neck+length+layer+2/3&source=bl&ots=YRM4VvBYZf&sig=ACfU3U13DlsYsJOp3sSzhrh3vI2LuarKow&hl=en&sa=X&ved=2ahUKEwjt25S-j-vgAhVtwcQBHcRwCFIQ6AEwB3oECAgQAQ#v=onepage&q=spine%20neck%20length%20layer%202%2F3&f=false (Yuste 2010 Dendritic spines, p. 149: 0.2-1.0um (or 0.1-2.0um))
neckDiam = 0.1 #(Yuste 2010 Dendritic spines, p.32: 0.1-0.5um)
headLen = 0.85 #These values of headlen and headdiam make up 0.5um3
headDiam = 0.85 #These values of headlen and headdiam make up 0.5um3
dodraw = 1
dendtree = 'apic'
spinelocations = '100-200'
if len(sys.argv) > 2:
  mutID = int(float(sys.argv[2]))
if len(sys.argv) > 3:
  stimfreq = float(sys.argv[3])
if len(sys.argv) > 4:
  stimN = int(float(sys.argv[4]))
if len(sys.argv) > 5:
  Nsyn = int(sys.argv[5])
if len(sys.argv) > 6:
  Ninputs = int(sys.argv[6])
if len(sys.argv) > 7:
  dendtree = sys.argv[7]
if len(sys.argv) > 8:
  spinelocations = sys.argv[8]
if len(sys.argv) > 9:
  rateE = float(sys.argv[9])
if len(sys.argv) > 10:
  Npulses = int(float(sys.argv[10]))
if len(sys.argv) > 11:
  ISIpairing = float(sys.argv[11])
if len(sys.argv) > 12:
  dtpulses = float(sys.argv[12])
if len(sys.argv) > 13:
  pulseamp = float(sys.argv[13])
if len(sys.argv) > 14:
  neckLen = float(sys.argv[14])
if len(sys.argv) > 15:
  neckDiam = float(sys.argv[15])
if len(sys.argv) > 16:
  Econ = float(sys.argv[16])
if len(sys.argv) > 17:
  weight_factor_NMDA = float(sys.argv[17])
if len(sys.argv) > 18:
  gNap = float(sys.argv[18])
if len(sys.argv) > 19:
  rdSeed = int(float(sys.argv[19]))
if len(sys.argv) > 20:
  dodraw = int(float(sys.argv[20]))

spineLocs = [int(float(x)) for x in spinelocations.split('-')]
if dendtree == 'apic':
  spineDendTree = 0
elif dendtree == 'dend' or dendtree == 'basal':
  spineDendTree = 1
else:
  print('dendtree '+dendtree+' not recognized!')
  sys.exit()

cellnames = ['cADpyr229_L23_PC_5ecbf9b163', 'cADpyr229_L23_PC_8ef1aa6602', 'cADpyr229_L23_PC_863902f300', 'cADpyr229_L23_PC_c292d67a2e', 'cADpyr229_L23_PC_c2e79db05a']
h("""
{load_file("stdlib.hoc")}
{load_file("stdrun.hoc")}
{load_file("L23_PC_cADpyr229_"""+str(icell+1)+"""/constants.hoc")}
{load_file("import3d.hoc")}
{load_file("morphology.hoc")}
{load_file("biophysics.hoc")}
{load_file("template.hoc")}
objref cvode
cvode = new CVode()
{cvode.active(1)}
{cvode.atol(0.00005)}

objref cell
objref time, voltage, iNMDAs_stim
cell = new """+cellnames[icell]+"""(0)
{voltage = new Vector()}
{time = new Vector()}
{iNMDAs_stim = new List()}
rdSeed = """+str(rdSeed)+"""
rateE = """+str(rateE)+"""
Napical = 0
Nbasal = 0
forsec cell.apical Napical = Napical + 1
forsec cell.basal Nbasal = Nbasal + 1
""")

areas = []
secnames = []
Nsegs = []
treenames = ['apic','dend']
treeIndStarts = [0, int(h.Napical)]
for i in range(0,int(h.Napical)+int(h.Nbasal)):
  if i < int(h.Napical):
    secname = 'apic['+str(i)+']'
  else:
    secname = 'dend['+str(i-int(h.Napical))+']'
  h("""myarea = 0
cell."""+secname+""" myNseg = nseg
for (i=0; i<myNseg; i+=1) { cell."""+secname+""" myarea = myarea + area((0.5+i)/myNseg) }
""")
  areas.append(h.myarea)
  secnames.append(secname)
  Nsegs.append(int(h.myNseg))
area_apic = sum(areas[0:int(h.Napical)])
area_dend = sum(areas[int(h.Napical):])
ps = [x/(area_apic+area_dend) for x in areas]
cumps = [sum(ps[0:1+i]) for i in range(0,len(ps))]

syncompsStim = []
synsecsStim = []
synxsStim = []
synxsegsStim = []
NspinesSet = 0
isyn = 0
while NspinesSet < Nsyn:
  if isyn%1000000 == 999999:
    print("Tried a million times...")
  isyn = isyn + 1
  r = rand()
  xseg = rand()
  isec = next((i for i in range(0,len(ps)) if cumps[i] > r))
  isecintree = isec if spineDendTree == 0 else isec - int(h.Napical)
  itree = int(isec >= int(h.Napical))
  if itree != spineDendTree:
    continue
  h('access cell.'+treenames[spineDendTree]+'['+str(isecintree)+']')
  mydist = h.distance(xseg)
  if mydist < spineLocs[0] or mydist > spineLocs[1]:
    continue
  syncompsStim.append(spineDendTree)
  synsecsStim.append(isecintree)
  synxsStim.append(xseg)
  synxsegsStim.append(int(synxsStim[-1]*Nsegs[isec]))
  NspinesSet = NspinesSet + 1
print("Spine locations determined after "+str(isyn)+" trials")

h("""
objref synlist, preconlist, nilstim, somastimlist
synlist = new List()
preconlist = new List()
somastimlist = new List()
""")

tstart_samp = 4500

for istim in range(0,stimN):
  for ipulse in range(0,Npulses):
    h('cell.soma somastimlist.append(new IClamp(0.5))')
    h('somastimlist.o(somastimlist.count()-1).dur = '+str(pulsedur))
    h('somastimlist.o(somastimlist.count()-1).amp = '+str(pulseamp))
    h('somastimlist.o(somastimlist.count()-1).del = '+str(tstart_samp+istim*(1000.0/stimfreq)+ISIpairing+ipulse*dtpulses))

#Poisson + stimulus inputs, exc.
h("SynList_stimulusSynstart = synlist.count()")
h("create spineNeck["+str(Nsyn)+"]")
h("create spineHead["+str(Nsyn)+"]")
stimcomps = []
stimdistances = []
stimsecs = []
stimxs = []
for isyn in range(0,Nsyn):
  h('spineNeck['+str(isyn)+'].L = '+str(neckLen))
  h('spineNeck['+str(isyn)+'].diam = '+str(neckLen))
  h('spineNeck['+str(isyn)+'].nseg = 1')
  h('spineNeck['+str(isyn)+'].cm = cell.soma.cm')
  h('spineNeck['+str(isyn)+'].Ra = cell.soma.Ra')
  h('spineHead['+str(isyn)+'].L = '+str(headLen))
  h('spineHead['+str(isyn)+'].diam = '+str(headLen))
  h('spineHead['+str(isyn)+'].nseg = 1')
  h('spineHead['+str(isyn)+'].cm = cell.soma.cm')
  h('spineHead['+str(isyn)+'].Ra = cell.soma.Ra')
  h('spineHead['+str(isyn)+'] insert Nap_Et2')
  h('spineHead['+str(isyn)+'] gNap_Et2bar_Nap_Et2 = '+str(gNap))
  icomp = synsecsStim[isyn]
  isec = icomp if icomp < int(h.Napical) else icomp - int(h.Napical)
  secname = 'apic' if icomp < int(h.Napical) else 'dend'
  x = synxsStim[isyn]
  h("""
cell."""+secname+"""["""+str(isec)+"""] connect spineNeck["""+str(isyn)+"""](0), """+str(x)+"""
spineNeck["""+str(isyn)+"""] connect spineHead["""+str(isyn)+"""](0), 1
access spineHead["""+str(isyn)+"""]
""")
  for iinput in range(0,Ninputs):
    h("""
{spineHead["""+str(isyn)+"""](0) synlist.append(new ProbAMPANMDA_EMST(0.5))}
iloc = synlist.count()-1
synlist.o[iloc].gmax = """+str(Econ)+"""
synlist.o[iloc].tau_r_AMPA = 0.3
synlist.o[iloc].tau_d_AMPA = 3
synlist.o[iloc].tau_r_NMDA = 2
synlist.o[iloc].tau_d_NMDA = 65
synlist.o[iloc].e = 0
synlist.o[iloc].Dep = 670
synlist.o[iloc].Use = 0.5
synlist.o[iloc].Fac = 17
synlist.o[iloc].u0 = 0
synlist.o[iloc].weight_factor_NMDA = """+str(weight_factor_NMDA)+"""
{preconlist.append(new NetCon(nilstim, synlist.o[iloc]))}
//print "len(preconlist) = ", preconlist.count(), ", len(synlist) = ", synlist.count()
preconi = preconlist.count()-1 //connection index
preconlist.o[preconi].weight = 1.0
preconlist.o[preconi].delay = 0
""")
  stimcomps.append(int(icomp >= int(h.Napical)))
  stimdistances.append(h.distance(x,sec=h.cell.apic[isec] if icomp < int(h.Napical) else h.cell.dend[isec]))
  stimsecs.append(isec)
  stimxs.append(x)

#LOAD MUTATIONS
mutArray = ['gCa_HVAbar', 'gCa_LVAstbar', 'gIhbar', 'gK_Pstbar', 'gSK_E2bar',  'g', 'gNaTs2_tbar','gImbar','gK_Tstbar','gSKv3_1bar','decay', 'gamma', 'gNap_Et2bar']
mutSuffs = ['Ca_HVA', 'Ca_LVAst', 'Ih', 'K_Pst', 'SK_E2', 'pas', 'NaTs2_t', 'Im', 'K_Tst', 'SKv3_1', 'CaDynamics_E2', 'CaDynamics_E2', 'Nap_Et2']
mutCoeffs = [0.75, 0.8, 0.9, 1.1, 1.2, 1.25]
mutText = ''
if mutID > 0:
  mutVar = mutArray[int((mutID-1)/6)]
  mutSuff = mutSuffs[int((mutID-1)/6)]
  mutCoeff = mutCoeffs[int((mutID-1)%6)]
  mutText = """forall if(ismembrane(\""""+mutSuff+"""\")) { """+mutVar+"""_"""+mutSuff+""" = """+str(mutCoeff)+"""*"""+mutVar+"""_"""+mutSuff+""" }"""
  print(mutText)
  h(mutText)
  if mutText.find('NaTs2_t') > -1:
    mutText2 = mutText.replace('NaTs2_t','NaTa_t')
    print(mutText2)
    h(mutText2)

h("""
access cell.soma
//time.record(&t, 0.1)
//{voltage.record(&v(0.5), 0.1)}
{cell.soma cvode.record(&v(0.5),voltage,time)}
""")
for isyn in range(0,Nsyn):
  icomp = synsecsStim[isyn]
  isec = icomp if icomp < int(h.Napical) else icomp - int(h.Napical)
  secname = 'apic' if icomp < int(h.Napical) else 'dend'
  x = synxsStim[isyn]
  for iinput in range(0,Ninputs):
    h('{iNMDAs_stim.append(new Vector())}')
    #h('{cell.'+secname+'['+str(isec)+'] cvode.record(&synlist.o[SynList_stimulusSynstart+'+str(isyn)+'].i_NMDA,iNMDAs_stim.o['+str(isyn)+'],time)}')
    h('{spineHead['+str(isyn)+'] cvode.record(&synlist.o[SynList_stimulusSynstart+'+str(isyn*Ninputs + iinput)+'].i_NMDA,iNMDAs_stim.o['+str(isyn*Ninputs + iinput)+'],time)}')

h("""
{objref fih,preTrainList,rds1}
{preTrainList = new List()}
{rds1 = new Random(1000*rdSeed+i)}//random for presynaptic trains
if (rateE > 0) {rds1.negexp(1/rateE)} else {rds1.negexp(1e8)}
""")

h("""
proc myqueue() {local isyn
  //Exc., Poisson + Stimulated
  for (isyn=0; isyn<"""+str(Nsyn*Ninputs)+"""; isyn+=1) {
    {preTrainList.append(new Vector())}
    if (rateE > 0) {
      pst=0 //presynaptic spike time
      while(pst < tstop){
        //print "{preconlist.o[", SynList_stimulusSynstart+isyn,"].event(pst), len(preconlist) = ", preconlist.count(), ", Nsyn = ", """+str(Nsyn)+""", ", len(synsecsStim) = """+str(len(synsecsStim))+""" "
        pst+= 1000*rds1.repick()
        {preTrainList.o[preTrainList.count()-1].append(pst)}
        {preconlist.o[SynList_stimulusSynstart+isyn].event(pst)}
      }
    }
    for (istim=0; istim<"""+str(stimN)+"""; istim+=1) {
      //tnow = """+str(tstart_samp)+"""+istim*"""+str(1000.0/stimfreq)+"""
      //print "event tnow = ", tnow
      {preTrainList.o[preTrainList.count()-1].append("""+str(tstart_samp)+"""+istim*"""+str(1000.0/stimfreq)+""")}
      {preconlist.o[SynList_stimulusSynstart+isyn].event("""+str(tstart_samp)+"""+istim*"""+str(1000.0/stimfreq)+""")}
    }
  }
}
""")

print("tstop = "+str(tstop))
timenow = time.time()
h("""
{fih = new FInitializeHandler("myqueue()")}
tstop = """+str(tstop)+"""
v0 = -75
init()
print "Starting simulation..."
run()""")
print('Simulation done in '+str(time.time()-timenow)+' seconds')

dt_samp = 1.0
Nsamp = 15000
#  Npulses = int(float(sys.argv[6]))
#  ISIpairing = float(sys.argv[7])
#  dtpulses = float(sys.argv[8])
scipy.io.savemat("noisy_icell"+str(icell)+"_imutc"+str(mutID)+"_n"+str(stimN)+"_"+str(stimfreq)+"_neckLen"+str(neckLen)+"_neckDiam"+str(neckDiam)+"_Nsyn"+str(Nsyn)+"_Ninputs"+str(Ninputs)+dendtree+spinelocations+"_Econ"+str(Econ)+'_wNMDA'+str(weight_factor_NMDA)+"_gNap"+str(gNap)+"_rateE"+str(rateE)+"_Npulses"+str(Npulses)+"_ISI"+str(ISIpairing)+"_dtpulses"+str(dtpulses)+"_pulseamp"+str(pulseamp)+"_seed"+str(rdSeed)+".mat",
                 {'tstart_samp': tstart_samp, 'dt_samp': dt_samp, 'Nsamp': Nsamp, 
                  'DATA': [mytools.interpolate(array(h.time),array(h.iNMDAs_stim[isyn]),[tstart_samp+i*dt_samp for i in range(0,Nsamp)]) for isyn in range(0,len(h.iNMDAs_stim))],
                  'stimdistances': stimdistances, 'stimcomps': stimcomps, 'stimsecs': stimsecs, 'stimxs': stimxs, 'vsoma': array(h.voltage), 'times': array(h.time), 'mutText': mutText})

if dodraw:
  f,axarr = subplots(1,1)
  axarr.plot(array(h.time),array(h.voltage))
  f.savefig("vsoma_icell"+str(icell)+"_imutc"+str(mutID)+"_n"+str(stimN)+"_"+str(stimfreq)+"_neckLen"+str(neckLen)+"_neckDiam"+str(neckDiam)+"_Nsyn"+str(Nsyn)+"_Ninputs"+str(Ninputs)+dendtree+spinelocations+"_Econ"+str(Econ)+'_wNMDA'+str(weight_factor_NMDA)+"_gNap"+str(gNap)+"_rateE"+str(rateE)+"_Npulses"+str(Npulses)+"_ISI"+str(ISIpairing)+"_dtpulses"+str(dtpulses)+"_pulseamp"+str(pulseamp)+"_seed"+str(rdSeed)+".eps")


