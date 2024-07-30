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
#mutArray = ['gCa_HVAbar', 'gCa_LVAstbar', 'gIhbar', 'gK_Pstbar', 'gSK_E2bar',  'g', 'gNaTs2_tbar','gImbar','gK_Tstbar','gSKv3_1bar','decay', 'gamma', 'gNap_Et2bar']
#Tuomo: Removed the fifth last one ('gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gImbar_Im,gK_Pstbar_K_Pst', [1.116,1.112,1.039,1.102]) to be consistent number with the mutcombs in spineNg5d)
mutSuffs = ['Ca_HVA', 'Ca_LVAst', 'Ih', 'K_Pst', 'SK_E2', 'pas', 'NaTs2_t', 'Im', 'K_Tst', 'SKv3_1', 'CaDynamics_E2', 'CaDynamics_E2', 'Nap_Et2']
#mutArray = ['gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gImbar_Im,gK_Pstbar_K_Pst', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gImbar_Im,gK_Pstbar_K_Pst', 'gIhbar_Ih,gNaTa_t_tbar_NaTa_t,gNaTs2_t_tbar_NaTs2_t,gNap_Et2bar_Nap_Et2,g_pas', 'gCa_HVAbar_Ca_HVA,gIhbar_Ih,gNap_Et2bar_Nap_Et2', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst,gNaTa_t_tbar_NaTa_t,gNaTs2_t_tbar_NaTs2_t,gNap_Et2bar_Nap_Et2,g_pas', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst,gNap_Et2bar_Nap_Et2', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gImbar_Im,gK_Pstbar_K_Pst', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gImbar_Im,gK_Pstbar_K_Pst', 'gIhbar_Ih,gNaTa_t_tbar_NaTa_t,gNaTs2_t_tbar_NaTs2_t,gNap_Et2bar_Nap_Et2,g_pas', 'gCa_HVAbar_Ca_HVA,gIhbar_Ih,gNap_Et2bar_Nap_Et2', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst,gNaTa_t_tbar_NaTa_t,gNaTs2_t_tbar_NaTs2_t,gNap_Et2bar_Nap_Et2,g_pas', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst,gNap_Et2bar_Nap_Et2', 'gCa_HVAbar_Ca_HVA,gImbar_Im,gK_Pstbar_K_Pst', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gImbar_Im,gK_Pstbar_K_Pst', 'gIhbar_Ih,gNaTa_t_tbar_NaTa_t,gNaTs2_t_tbar_NaTs2_t,gNap_Et2bar_Nap_Et2,g_pas', 'gCa_HVAbar_Ca_HVA,gIhbar_Ih,gNap_Et2bar_Nap_Et2', 'gCa_HVAbar_Ca_HVA,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst,gNaTa_t_tbar_NaTa_t,gNaTs2_t_tbar_NaTs2_t,gNap_Et2bar_Nap_Et2,g_pas', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst,gNap_Et2bar_Nap_Et2', 'gCa_HVAbar_Ca_HVA,gImbar_Im,gK_Pstbar_K_Pst', 'gIhbar_Ih,gNaTa_t_tbar_NaTa_t,gNaTs2_t_tbar_NaTs2_t,gNap_Et2bar_Nap_Et2,g_pas', 'gCa_HVAbar_Ca_HVA,gIhbar_Ih,gNap_Et2bar_Nap_Et2', 'gCa_HVAbar_Ca_HVA,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst,gNaTa_t_tbar_NaTa_t,gNaTs2_t_tbar_NaTs2_t,gNap_Et2bar_Nap_Et2,g_pas', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst,gNap_Et2bar_Nap_Et2'] #old one, with zero Ih coeffs due to wrong numbering of igenes in the old extract_params_from_genes_all_review.py
#mutArray = ['gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst', 'gIhbar_Ih,gNaTa_t_tbar_NaTa_t,gNaTs2_t_tbar_NaTs2_t,gNap_Et2bar_Nap_Et2,g_pas', 'gCa_HVAbar_Ca_HVA,gK_Pstbar_K_Pst,gNap_Et2bar_Nap_Et2', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst,gNaTa_t_tbar_NaTa_t,gNaTs2_t_tbar_NaTs2_t,gNap_Et2bar_Nap_Et2,g_pas', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst,gNap_Et2bar_Nap_Et2', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gImbar_Im,gK_Pstbar_K_Pst', 'gIhbar_Ih,gNaTa_t_tbar_NaTa_t,gNaTs2_t_tbar_NaTs2_t,gNap_Et2bar_Nap_Et2,g_pas', 'gCa_HVAbar_Ca_HVA,gK_Pstbar_K_Pst,gNap_Et2bar_Nap_Et2', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst,gNaTa_t_tbar_NaTa_t,gNaTs2_t_tbar_NaTs2_t,gNap_Et2bar_Nap_Et2,g_pas', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gImbar_Im,gK_Pstbar_K_Pst,gNap_Et2bar_Nap_Et2', 'gCa_HVAbar_Ca_HVA,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gImbar_Im,gK_Pstbar_K_Pst', 'gIhbar_Ih,gNaTa_t_tbar_NaTa_t,gNaTs2_t_tbar_NaTs2_t,gNap_Et2bar_Nap_Et2,g_pas', 'gCa_HVAbar_Ca_HVA,gK_Pstbar_K_Pst,gNap_Et2bar_Nap_Et2', 'gCa_HVAbar_Ca_HVA,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst,gNaTa_t_tbar_NaTa_t,gNaTs2_t_tbar_NaTs2_t,gNap_Et2bar_Nap_Et2,g_pas', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gImbar_Im,gK_Pstbar_K_Pst,gNap_Et2bar_Nap_Et2', 'gCa_HVAbar_Ca_HVA,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst', 'gIhbar_Ih,gNaTa_t_tbar_NaTa_t,gNaTs2_t_tbar_NaTs2_t,gNap_Et2bar_Nap_Et2,g_pas', 'gCa_HVAbar_Ca_HVA,gK_Pstbar_K_Pst,gNap_Et2bar_Nap_Et2', 'gCa_HVAbar_Ca_HVA,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst,gNaTa_t_tbar_NaTa_t,gNaTs2_t_tbar_NaTs2_t,gNap_Et2bar_Nap_Et2,g_pas', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gImbar_Im,gK_Pstbar_K_Pst,gNap_Et2bar_Nap_Et2'] #wrong Na current assigned for SCN9A, wrong name of Na currents
mutArray = ['gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst', 'gIhbar_Ih,gNaTa_tbar_NaTa_t,gNaTs2_tbar_NaTs2_t,g_pas', 'gCa_HVAbar_Ca_HVA,gK_Pstbar_K_Pst,gNaTa_tbar_NaTa_t,gNaTs2_tbar_NaTs2_t', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst,gNaTa_tbar_NaTa_t,gNaTs2_tbar_NaTs2_t,g_pas', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst,gNaTa_tbar_NaTa_t,gNaTs2_tbar_NaTs2_t', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gImbar_Im,gK_Pstbar_K_Pst', 'gIhbar_Ih,gNaTa_tbar_NaTa_t,gNaTs2_tbar_NaTs2_t,g_pas', 'gCa_HVAbar_Ca_HVA,gK_Pstbar_K_Pst,gNaTa_tbar_NaTa_t,gNaTs2_tbar_NaTs2_t', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst,gNaTa_tbar_NaTa_t,gNaTs2_tbar_NaTs2_t,g_pas', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gImbar_Im,gK_Pstbar_K_Pst,gNaTa_tbar_NaTa_t,gNaTs2_tbar_NaTs2_t', 'gCa_HVAbar_Ca_HVA,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gImbar_Im,gK_Pstbar_K_Pst', 'gIhbar_Ih,gNaTa_tbar_NaTa_t,gNaTs2_tbar_NaTs2_t,g_pas', 'gCa_HVAbar_Ca_HVA,gK_Pstbar_K_Pst,gNaTa_tbar_NaTa_t,gNaTs2_tbar_NaTs2_t', 'gCa_HVAbar_Ca_HVA,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst,gNaTa_tbar_NaTa_t,gNaTs2_tbar_NaTs2_t,g_pas', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gImbar_Im,gK_Pstbar_K_Pst,gNaTa_tbar_NaTa_t,gNaTs2_tbar_NaTs2_t', 'gCa_HVAbar_Ca_HVA,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst', 'gIhbar_Ih,gNaTa_tbar_NaTa_t,gNaTs2_tbar_NaTs2_t,g_pas', 'gCa_HVAbar_Ca_HVA,gK_Pstbar_K_Pst,gNaTa_tbar_NaTa_t,gNaTs2_tbar_NaTs2_t', 'gCa_HVAbar_Ca_HVA,gIhbar_Ih,gImbar_Im,gK_Pstbar_K_Pst,gNaTa_tbar_NaTa_t,gNaTs2_tbar_NaTs2_t,g_pas', 'gCa_HVAbar_Ca_HVA,gCa_LVAstbar_Ca_LVAst,gImbar_Im,gK_Pstbar_K_Pst,gNaTa_tbar_NaTa_t,gNaTs2_tbar_NaTs2_t']
#mutCoeffs = [[1.068,1.023,1.077,1.087], [1.126,1.138,1.077,1.125], [1.136,0.829,0.829,1.368,0.841], [1.126,0.000,1.240], [1.068,1.023,1.136,1.077,1.087,0.829,0.829,1.368,0.841], [1.126,1.138,0.000,1.077,1.125,1.240], [1.027,1.029,1.037,1.076], [1.116,1.112,1.039,1.102], [1.061,0.915,0.915,1.275,0.897], [1.116,0.000,1.174], [1.027,1.029,1.061,1.037,1.076,0.915,0.915,1.275,0.897], [1.116,1.112,0.000,1.039,1.102,1.174], [1.074,1.077,1.087], [1.126,1.138,1.077,1.125], [1.136,0.829,0.829,1.368,0.841], [1.126,0.000,1.240], [1.074,1.136,1.077,1.087,0.829,0.829,1.368,0.841], [1.126,1.138,0.000,1.077,1.125,1.240], [1.017,1.037,1.076], [1.061,0.915,0.915,1.275,0.897], [1.116,0.000,1.174], [1.017,1.061,1.037,1.076,0.915,0.915,1.275,0.897], [1.116,1.112,0.000,1.039,1.102,1.174]] #old one, with zero Ih coeffs due to wrong numbering of igenes in the old extract_params_from_genes_all_review.py
#mutCoeffs = [[1.068,1.023,1.136,1.077,1.087], [1.126,1.138,1.056,1.077,1.125], [1.136,0.829,0.829,1.368,0.841], [1.126,1.080,1.240], [1.068,1.023,1.136,1.077,1.087,0.829,0.829,1.368,0.841], [1.126,1.138,1.056,1.077,1.102,1.240], [1.027,1.029,1.061,1.037,1.076], [1.116,1.112,1.039,1.102], [1.061,0.915,0.915,1.275,0.897], [1.116,1.059,1.174], [1.027,1.029,1.061,1.037,1.076,0.915,0.915,1.275,0.897], [1.116,1.112,1.039,1.081,1.174], [1.074,1.136,1.077,1.087], [1.126,1.138,1.077,1.125], [1.136,0.829,0.829,1.368,0.841], [1.126,1.080,1.240], [1.074,1.136,1.077,1.087,0.829,0.829,1.368,0.841], [1.126,1.138,1.077,1.102,1.240], [1.017,1.061,1.037,1.076], [1.061,0.915,0.915,1.275,0.897], [1.116,1.059,1.174], [1.017,1.061,1.037,1.076,0.915,0.915,1.275,0.897], [1.116,1.112,1.039,1.081,1.174]] #wrong Na current assigned for SCN9A, wrong name of Na currents
mutCoeffs = [[1.068,1.023,1.136,1.077,1.087], [1.126,1.138,1.056,1.077,1.125], [1.136,1.099,1.099,0.841], [1.126,1.080,1.240,1.240], [1.068,1.023,1.136,1.077,1.087,1.099,1.099,0.841], [1.126,1.138,1.056,1.077,1.102,1.240,1.240], [1.027,1.029,1.061,1.037,1.076], [1.116,1.112,1.039,1.102], [1.061,1.095,1.095,0.897], [1.116,1.059,1.174,1.174], [1.027,1.029,1.061,1.037,1.076,1.095,1.095,0.897], [1.116,1.112,1.039,1.081,1.174,1.174], [1.074,1.136,1.077,1.087], [1.126,1.138,1.077,1.125], [1.136,1.099,1.099,0.841], [1.126,1.080,1.240,1.240], [1.074,1.136,1.077,1.087,1.099,1.099,0.841], [1.126,1.138,1.077,1.102,1.240,1.240], [1.017,1.061,1.037,1.076], [1.061,1.095,1.095,0.897], [1.116,1.059,1.174,1.174], [1.017,1.061,1.037,1.076,1.095,1.095,0.897], [1.116,1.112,1.039,1.081,1.174,1.174]]

mutText = ''
if mutID > 0:
  mutVars = mutArray[mutID]
  mutCoeffs = mutCoeffs[mutID]
  mutVarsSpl = mutVars.split(',')
  for imutvar in range(0,len(mutVarsSpl)):
    mutVar = mutVarsSpl[imutvar]
    mutCoeff = mutCoeffs[imutvar]
    suffsFound = [mutSuffs[i] for i in range(0,len(mutSuffs)) if mutSuffs[i] in mutVar]
    if len(suffsFound) == 0:
      print('Suffix not found, mutVar='+mutVar)
      continue
    mutSuff = suffsFound[0]
    isuffstart = mutVar.find(mutSuff)
    mutText = """forall if(ismembrane(\""""+mutSuff+"""\")) { """+mutVar+""" = """+str(mutCoeff)+"""*"""+mutVar+""" }"""
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
  scipy.io.savemat("somaticDC_icell"+str(icell)+"_atol"+str(atol)+"_imutCMcomb"+str(mutID)+"_stimAmp"+str(stimAmp)+".mat", {'vsoma': array(h.voltage), 'times': array(h.time), 'mutText': mutText})
else:
  scipy.io.savemat("somaticDC_icell"+str(icell)+"_imutCMcomb"+str(mutID)+"_stimAmp"+str(stimAmp)+".mat", {'vsoma': array(h.voltage), 'times': array(h.time), 'mutText': mutText})
