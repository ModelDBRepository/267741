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

mutID = 0
stimAmp = 0.0
stimOnset = 1000
stimDur = 19000
dodraw = 1
atol = 0.00005
if len(sys.argv) > 2:
  mutID = int(float(sys.argv[2]))
if len(sys.argv) > 3:
  stimAmp = float(sys.argv[3])
if len(sys.argv) > 4:
  atol = float(sys.argv[4])



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
{cvode.atol("""+str(atol)+""")}

objref cell
objref time, voltage, iNMDAs_stim
cell = new """+cellnames[icell]+"""(0)
{voltage = new Vector()}
{time = new Vector()}
{iNMDAs_stim = new List()}
Napical = 0
Nbasal = 0
forsec cell.apical Napical = Napical + 1
forsec cell.basal Nbasal = Nbasal + 1
""")

h('objref somastim')
h('cell.soma somastim = new IClamp(0.5)')
h('somastim.dur = '+str(stimDur))
h('somastim.amp = '+str(stimAmp))
h('somastim.del = '+str(stimOnset))


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

timenow = time.time()
h("""
tstop = """+str(tstop)+"""
v0 = -75
init()
print "Starting simulation..."
run()""")
print('Simulation done in '+str(time.time()-timenow)+' seconds')

if atol != 0.00005:
  scipy.io.savemat("somaticDC_icell"+str(icell)+"_atol"+str(atol)+"_imutc"+str(mutID)+"_stimAmp"+str(stimAmp)+".mat", {'vsoma': array(h.voltage), 'times': array(h.time), 'mutText': mutText})
else:
  scipy.io.savemat("somaticDC_icell"+str(icell)+"_imutc"+str(mutID)+"_stimAmp"+str(stimAmp)+".mat", {'vsoma': array(h.voltage), 'times': array(h.time), 'mutText': mutText})
