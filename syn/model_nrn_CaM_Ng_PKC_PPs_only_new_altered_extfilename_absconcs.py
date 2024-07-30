from neuron import h, rxd
from pylab import *
from matplotlib import pyplot
import scipy.io
import time
import re
import mytools
from os.path import exists

h.load_file('stdrun.hoc')

dend = h.Section(name='dend')
dend.L=1
dend.diam=0.79788
cyt = rxd.Region([dend], name='cyt', nrn_region='i')

mesh_input_file = open('mesh_general.out','r')
mesh_firstline = mesh_input_file.readline()
mesh_secondline = mesh_input_file.readline()
mesh_values = mesh_secondline.split()
my_volume = float(mesh_values[-2])*1e-15 #litres
mesh_input_file.close()

Duration = 1000
tolerance = 1e-9
initfile = ''
alteredk = []
alteredks = []
altered_factor = 1.0
tolstochange = 'Gs,Calbin,CalbinC,PMCA,PMCACa,PP2B,Gi,GiaGDP,Gibg,Gsbg,GsaGDP,CaOut,CaOutLeak,Leak,PDE1,Ng,NgCaM,NgpCaM,CaM,CaMCa2,CaMCa4,fixedbuffer,fixedbufferCa,DGL,CaDGL'.split(',')
tolschange = [float(x) for x in '3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,2,2,1,1'.split(',')]

if len(sys.argv) > 1:
  Duration = int(sys.argv[1])
if len(sys.argv) > 2:
  tolerance = float(sys.argv[2])
if len(sys.argv) > 3:
  speciesgiven_str = sys.argv[3]
  speciesgiven = speciesgiven_str.split(',')
if len(sys.argv) > 4:
  concsgiven_str = sys.argv[4]
  concsgiven = [float(x) for x in concsgiven_str.split(',')]
if len(sys.argv) > 5:
  alteredk = sys.argv[5]
  alteredks = [int(x) for x in alteredk.split(',')]
if len(sys.argv) > 6:
  alteredk_factor = sys.argv[6]
  alteredk_factors = [float(x) for x in alteredk_factor.split(',')]
filename = 'nrn_tstop'+str(Duration)+'_tol'+str(tolerance)+'.mat'
if len(sys.argv) > 7:
  filename = sys.argv[7]
if len(sys.argv) > 8:
  toltochange = 'Gs,Calbin,CalbinC,PMCA,PMCACa,PP2B,Gi,GiaGDP,Gibg,Gsbg,GsaGDP,CaOut,CaOutLeak,Leak,PDE1,Ng,NgCaM,NgpCaM,CaM,CaMCa2,CaMCa4,fixedbuffer,fixedbufferCa,DGL,CaDGL,'+sys.argv[8]
  tolstochange = toltochange.split(',')
if len(sys.argv) > 9:
  tolchange = '3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,2,2,1,1,'+sys.argv[9]
  tolschange = [float(x) for x in tolchange.split(',')]

initvalues = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0]

species = ['Ca', 'Ng', 'NgCaM', 'Ngp', 'NgpCaM', 'NgPKCt', 'NgCaMPKCt', 'NgPKCp', 'NgCaMPKCp', 'NgCaMCa2', 'NgpCaMCa2', 'NgCaMCa2PKCt', 'NgCaMCa2PKCp', 'NgCaMCa3', 'NgpCaMCa3', 'NgCaMCa3PKCt', 'NgCaMCa3PKCp', 'NgCaMCa4', 'NgpCaMCa4', 'NgCaMCa4PKCt', 'NgCaMCa4PKCp', 'NgpPP1', 'NgpPP2A', 'NgpPP2BCaMCa4', 'PP1', 'PP2A', 'CaM', 'CaMCa2', 'CaMCa3', 'CaMCa4', 'PP2B', 'PP2BCaM', 'PP2BCaMCa2', 'PP2BCaMCa3', 'PP2BCaMCa4', 'PKC', 'PKCCa', 'PKCt', 'PKCp']

print("my_volume = "+str(my_volume)+" l ?= "+str(dend.L*(dend.diam/2)**2*3.14159265358)+" um3")
if len(initfile) > 3 and initfile != 'None':
  DATA_init = scipy.io.loadmat(initfile)
  for mykey in DATA_init.keys():
    if mykey[0:2] != '__':
      DATA_init[mykey] = DATA_init[mykey][0]
  for ispec in range(0,len(species)):
    initvalues[ispec] = DATA_init[species][-1]
tolscales = [1.0 for i in range(0,len(species))]
for ispec in range(0,len(species)):
  for ispecgiven in range(0,len(speciesgiven)):
    if re.match(speciesgiven[ispecgiven]+'$',species[ispec]):
      initvalues[ispec] = concsgiven[ispecgiven]
  for itol in range(0,len(tolstochange)):
    if tolstochange[itol] == species[ispec]:
      tolscales[ispec] = tolscales[ispec]*10**(-tolschange[itol])
specs = []

for ispec in range(0,len(species)):
  specs.append(rxd.Species(cyt, name='spec'+str(ispec), charge=0, initial=initvalues[ispec], atolscale=tolscales[ispec]))
  if tolscales[ispec] != 1.0:
    print('spec'+str(ispec)+' ('+species[ispec]+') atolscale = '+str(tolscales[ispec]))
ks = [1.0]*94
ks[0]   = 17000.0 # CaM + Ca*2 <-> CaMCa2 (forward)
ks[1]   = 0.035   # CaM + Ca*2 <-> CaMCa2 (backward)
ks[2]   = 14.0    # CaMCa2 + Ca <-> CaMCa3 (forward)
ks[3]   = 0.228   # CaMCa2 + Ca <-> CaMCa3 (backward)
ks[4]   = 26.0    # CaMCa3 + Ca <-> CaMCa4 (forward)
ks[5]   = 0.064   # CaMCa3 + Ca <-> CaMCa4 (backward)
ks[6]   = 0.46    # CaM + PP2B <-> PP2BCaM (forward)
ks[7]   = 1.2e-06 # CaM + PP2B <-> PP2BCaM (backward)
ks[8]   = 0.46    # CaMCa2 + PP2B <-> PP2BCaMCa2 (forward)
ks[9]   = 1.2e-06 # CaMCa2 + PP2B <-> PP2BCaMCa2 (backward)
ks[10]  = 4.6     # CaMCa3 + PP2B <-> PP2BCaMCa3 (forward)
ks[11]  = 1.2e-06 # CaMCa3 + PP2B <-> PP2BCaMCa3 (backward)
ks[12]  = 46.0    # CaMCa4 + PP2B <-> PP2BCaMCa4 (forward)
ks[13]  = 1.2e-06 # CaMCa4 + PP2B <-> PP2BCaMCa4 (backward)
ks[14]  = 17000.0 # PP2BCaM + Ca*2 <-> PP2BCaMCa2 (forward)
ks[15]  = 0.035   # PP2BCaM + Ca*2 <-> PP2BCaMCa2 (backward)
ks[16]  = 14.0    # PP2BCaMCa2 + Ca <-> PP2BCaMCa3 (forward)
ks[17]  = 0.0228  # PP2BCaMCa2 + Ca <-> PP2BCaMCa3 (backward)
ks[18]  = 26.0    # PP2BCaMCa3 + Ca <-> PP2BCaMCa4 (forward)
ks[19]  = 0.0064  # PP2BCaMCa3 + Ca <-> PP2BCaMCa4 (backward)
ks[20]  = 15806.7 # NgCaM + Ca*2 <-> NgCaMCa2 (forward)
ks[21]  = 1.72117 # NgCaM + Ca*2 <-> NgCaMCa2 (backward)
ks[22]  = 14.0    # NgCaMCa2 + Ca <-> NgCaMCa3 (forward)
ks[23]  = 0.228   # NgCaMCa2 + Ca <-> NgCaMCa3 (backward)
ks[24]  = 26.0    # NgCaMCa3 + Ca <-> NgCaMCa4 (forward)
ks[25]  = 0.064   # NgCaMCa3 + Ca <-> NgCaMCa4 (backward)
ks[26]  = 15806.7 # NgpCaM + Ca*2 <-> NgpCaMCa2 (forward)
ks[27]  = 1.72117 # NgpCaM + Ca*2 <-> NgpCaMCa2 (backward)
ks[28]  = 14.0    # NgpCaMCa2 + Ca <-> NgpCaMCa3 (forward)
ks[29]  = 0.228   # NgpCaMCa2 + Ca <-> NgpCaMCa3 (backward)
ks[30]  = 26.0    # NgpCaMCa3 + Ca <-> NgpCaMCa4 (forward)
ks[31]  = 0.064   # NgpCaMCa3 + Ca <-> NgpCaMCa4 (backward)
ks[32]  = 28.0    # CaM + Ng <-> NgCaM (forward)
ks[33]  = 0.036   # CaM + Ng <-> NgCaM (backward)
ks[34]  = 0.0     # CaM + Ngp <-> NgpCaM (forward)
ks[35]  = 0.036   # CaM + Ngp <-> NgpCaM (backward)
ks[36]  = 2.0     # CaMCa2 + Ng <-> NgCaMCa2 (forward)
ks[37]  = 0.136   # CaMCa2 + Ng <-> NgCaMCa2 (backward)
ks[38]  = 0.0     # CaMCa2 + Ngp <-> NgpCaMCa2 (forward)
ks[39]  = 0.136   # CaMCa2 + Ngp <-> NgpCaMCa2 (backward)
ks[40]  = 2.0     # CaMCa3 + Ng <-> NgCaMCa3 (forward)
ks[41]  = 0.136   # CaMCa3 + Ng <-> NgCaMCa3 (backward)
ks[42]  = 0.0     # CaMCa3 + Ngp <-> NgpCaMCa3 (forward)
ks[43]  = 0.136   # CaMCa3 + Ngp <-> NgpCaMCa3 (backward)
ks[44]  = 2.0     # CaMCa4 + Ng <-> NgCaMCa4 (forward)
ks[45]  = 0.136   # CaMCa4 + Ng <-> NgCaMCa4 (backward)
ks[46]  = 0.0     # CaMCa4 + Ngp <-> NgpCaMCa4 (forward)
ks[47]  = 0.136   # CaMCa4 + Ngp <-> NgpCaMCa4 (backward)
ks[48]  = 13.3    # Ca + PKC <-> PKCCa (forward)
ks[49]  = 0.05    # Ca + PKC <-> PKCCa (backward)
ks[50]  = 0.0678  # Ng + PKCt <-> NgPKCt (forward)
ks[51]  = 0.002   # Ng + PKCt <-> NgPKCt (backward)
ks[52]  = 1.775   # NgPKCt --> Ngp + PKCt (forward)
ks[53]  = 0.0678  # NgCaM + PKCt <-> NgCaMPKCt (forward)
ks[54]  = 0.002   # NgCaM + PKCt <-> NgCaMPKCt (backward)
ks[55]  = 0.0678  # NgCaMCa2 + PKCt <-> NgCaMCa2PKCt (forward)
ks[56]  = 0.002   # NgCaMCa2 + PKCt <-> NgCaMCa2PKCt (backward)
ks[57]  = 0.0678  # NgCaMCa3 + PKCt <-> NgCaMCa3PKCt (forward)
ks[58]  = 0.002   # NgCaMCa3 + PKCt <-> NgCaMCa3PKCt (backward)
ks[59]  = 0.0678  # NgCaMCa4 + PKCt <-> NgCaMCa4PKCt (forward)
ks[60]  = 0.002   # NgCaMCa4 + PKCt <-> NgCaMCa4PKCt (backward)
ks[61]  = 1.775   # NgCaMPKCt --> NgpCaM + PKCt (forward)
ks[62]  = 1.775   # NgCaMCa2PKCt --> NgpCaMCa2 + PKCt (forward)
ks[63]  = 1.775   # NgCaMCa3PKCt --> NgpCaMCa3 + PKCt (forward)
ks[64]  = 1.775   # NgCaMCa4PKCt --> NgpCaMCa4 + PKCt (forward)
ks[65]  = 0.0678  # Ng + PKCp <-> NgPKCp (forward)
ks[66]  = 0.002   # Ng + PKCp <-> NgPKCp (backward)
ks[67]  = 1.775   # NgPKCp --> Ngp + PKCp (forward)
ks[68]  = 0.0678  # NgCaM + PKCp <-> NgCaMPKCp (forward)
ks[69]  = 0.002   # NgCaM + PKCp <-> NgCaMPKCp (backward)
ks[70]  = 0.0678  # NgCaMCa2 + PKCp <-> NgCaMCa2PKCp (forward)
ks[71]  = 0.002   # NgCaMCa2 + PKCp <-> NgCaMCa2PKCp (backward)
ks[72]  = 0.0678  # NgCaMCa3 + PKCp <-> NgCaMCa3PKCp (forward)
ks[73]  = 0.002   # NgCaMCa3 + PKCp <-> NgCaMCa3PKCp (backward)
ks[74]  = 0.0678  # NgCaMCa4 + PKCp <-> NgCaMCa4PKCp (forward)
ks[75]  = 0.002   # NgCaMCa4 + PKCp <-> NgCaMCa4PKCp (backward)
ks[76]  = 1.775   # NgCaMPKCp --> NgpCaM + PKCp (forward)
ks[77]  = 1.775   # NgCaMCa2PKCp --> NgpCaMCa2 + PKCp (forward)
ks[78]  = 1.775   # NgCaMCa3PKCp --> NgpCaMCa3 + PKCp (forward)
ks[79]  = 1.775   # NgCaMCa4PKCp --> NgpCaMCa4 + PKCp (forward)
ks[80]  = 2.5e-06 # Ngp --> Ng (forward)
ks[81]  = 2.5e-06 # NgpCaM --> NgCaM (forward)
ks[82]  = 2.5e-06 # NgpCaMCa2 --> NgCaMCa2 (forward)
ks[83]  = 2.5e-06 # NgpCaMCa3 --> NgCaMCa3 (forward)
ks[84]  = 2.5e-06 # NgpCaMCa4 --> NgCaMCa4 (forward)
ks[85]  = 0.875   # Ngp + PP1 <-> NgpPP1 (forward)
ks[86]  = 0.0014  # Ngp + PP1 <-> NgpPP1 (backward)
ks[87]  = 0.00035 # NgpPP1 --> Ng + PP1 (forward)
ks[88]  = 0.5     # Ngp + PP2A <-> NgpPP2A (forward)
ks[89]  = 0.005   # Ngp + PP2A <-> NgpPP2A (backward)
ks[90]  = 0.00015 # NgpPP2A --> Ng + PP2A (forward)
ks[91]  = 2.01    # Ngp + PP2BCaMCa4 <-> NgpPP2BCaMCa4 (forward)
ks[92]  = 0.008   # Ngp + PP2BCaMCa4 <-> NgpPP2BCaMCa4 (backward)
ks[93]  = 0.002   # NgpPP2BCaMCa4 --> Ng + PP2BCaMCa4 (forward)

for ialteredk in range(0,len(alteredks)):
  ks[alteredks[ialteredk]] = alteredk_factors[ialteredk]*ks[alteredks[ialteredk]]
reaction000 = rxd.Reaction(specs[26] + specs[0]*2, specs[27], ks[0], ks[1])
reaction001 = rxd.Reaction(specs[27] + specs[0], specs[28], ks[2], ks[3])
reaction002 = rxd.Reaction(specs[28] + specs[0], specs[29], ks[4], ks[5])
reaction003 = rxd.Reaction(specs[26] + specs[30], specs[31], ks[6], ks[7])
reaction004 = rxd.Reaction(specs[27] + specs[30], specs[32], ks[8], ks[9])
reaction005 = rxd.Reaction(specs[28] + specs[30], specs[33], ks[10], ks[11])
reaction006 = rxd.Reaction(specs[29] + specs[30], specs[34], ks[12], ks[13])
reaction007 = rxd.Reaction(specs[31] + specs[0]*2, specs[32], ks[14], ks[15])
reaction008 = rxd.Reaction(specs[32] + specs[0], specs[33], ks[16], ks[17])
reaction009 = rxd.Reaction(specs[33] + specs[0], specs[34], ks[18], ks[19])
reaction010 = rxd.Reaction(specs[2] + specs[0]*2, specs[9], ks[20], ks[21])
reaction011 = rxd.Reaction(specs[9] + specs[0], specs[13], ks[22], ks[23])
reaction012 = rxd.Reaction(specs[13] + specs[0], specs[17], ks[24], ks[25])
reaction013 = rxd.Reaction(specs[4] + specs[0]*2, specs[10], ks[26], ks[27])
reaction014 = rxd.Reaction(specs[10] + specs[0], specs[14], ks[28], ks[29])
reaction015 = rxd.Reaction(specs[14] + specs[0], specs[18], ks[30], ks[31])
reaction016 = rxd.Reaction(specs[26] + specs[1], specs[2], ks[32], ks[33])
reaction017 = rxd.Reaction(specs[26] + specs[3], specs[4], ks[34], ks[35])
reaction018 = rxd.Reaction(specs[27] + specs[1], specs[9], ks[36], ks[37])
reaction019 = rxd.Reaction(specs[27] + specs[3], specs[10], ks[38], ks[39])
reaction020 = rxd.Reaction(specs[28] + specs[1], specs[13], ks[40], ks[41])
reaction021 = rxd.Reaction(specs[28] + specs[3], specs[14], ks[42], ks[43])
reaction022 = rxd.Reaction(specs[29] + specs[1], specs[17], ks[44], ks[45])
reaction023 = rxd.Reaction(specs[29] + specs[3], specs[18], ks[46], ks[47])
reaction024 = rxd.Reaction(specs[0] + specs[35], specs[36], ks[48], ks[49])
reaction025 = rxd.Reaction(specs[1] + specs[37], specs[5], ks[50], ks[51])
reaction026 = rxd.Reaction(specs[5], specs[3] + specs[37], ks[52])
reaction027 = rxd.Reaction(specs[2] + specs[37], specs[6], ks[53], ks[54])
reaction028 = rxd.Reaction(specs[9] + specs[37], specs[11], ks[55], ks[56])
reaction029 = rxd.Reaction(specs[13] + specs[37], specs[15], ks[57], ks[58])
reaction030 = rxd.Reaction(specs[17] + specs[37], specs[19], ks[59], ks[60])
reaction031 = rxd.Reaction(specs[6], specs[4] + specs[37], ks[61])
reaction032 = rxd.Reaction(specs[11], specs[10] + specs[37], ks[62])
reaction033 = rxd.Reaction(specs[15], specs[14] + specs[37], ks[63])
reaction034 = rxd.Reaction(specs[19], specs[18] + specs[37], ks[64])
reaction035 = rxd.Reaction(specs[1] + specs[38], specs[7], ks[65], ks[66])
reaction036 = rxd.Reaction(specs[7], specs[3] + specs[38], ks[67])
reaction037 = rxd.Reaction(specs[2] + specs[38], specs[8], ks[68], ks[69])
reaction038 = rxd.Reaction(specs[9] + specs[38], specs[12], ks[70], ks[71])
reaction039 = rxd.Reaction(specs[13] + specs[38], specs[16], ks[72], ks[73])
reaction040 = rxd.Reaction(specs[17] + specs[38], specs[20], ks[74], ks[75])
reaction041 = rxd.Reaction(specs[8], specs[4] + specs[38], ks[76])
reaction042 = rxd.Reaction(specs[12], specs[10] + specs[38], ks[77])
reaction043 = rxd.Reaction(specs[16], specs[14] + specs[38], ks[78])
reaction044 = rxd.Reaction(specs[20], specs[18] + specs[38], ks[79])
reaction045 = rxd.Reaction(specs[3], specs[1], ks[80])
reaction046 = rxd.Reaction(specs[4], specs[2], ks[81])
reaction047 = rxd.Reaction(specs[10], specs[9], ks[82])
reaction048 = rxd.Reaction(specs[14], specs[13], ks[83])
reaction049 = rxd.Reaction(specs[18], specs[17], ks[84])
reaction050 = rxd.Reaction(specs[3] + specs[24], specs[21], ks[85], ks[86])
reaction051 = rxd.Reaction(specs[21], specs[1] + specs[24], ks[87])
reaction052 = rxd.Reaction(specs[3] + specs[25], specs[22], ks[88], ks[89])
reaction053 = rxd.Reaction(specs[22], specs[1] + specs[25], ks[90])
reaction054 = rxd.Reaction(specs[3] + specs[34], specs[23], ks[91], ks[92])
reaction055 = rxd.Reaction(specs[23], specs[1] + specs[34], ks[93])
vec_t = h.Vector()

vecs = []
vec_t = h.Vector()
vec_t.record(h._ref_t)
for ispec in range(0,len(species)):
  vecs.append(h.Vector())
  vecs[ispec].record(specs[ispec].nodes(dend)(0.5)[0]._ref_concentration)

cvode = h.CVode()
cvode.active(1)
hmax = cvode.maxstep(1000)
hmin = cvode.minstep(1e-10)
cvode.atol(tolerance)

h.finitialize(-65)
def set_param(param, val):
    param.nodes.value = val
    h.cvode.re_init()

timenow = time.time()
h.continuerun(Duration)
print("Simulation done in "+str(time.time()-timenow)+" seconds")
tvec = array(vec_t)
minDT = 1.0
lastt = -inf
itvec2 = []
for it in range(0,len(tvec)):
  if tvec[it] - lastt > minDT:
    itvec2.append(it)
    lastt = tvec[it]

headers = [ 'tvec', 'Ca', 'Ng', 'NgCaM', 'Ngp', 'NgpCaM', 'NgPKCt', 'NgCaMPKCt', 'NgPKCp', 'NgCaMPKCp', 'NgCaMCa2', 'NgpCaMCa2', 'NgCaMCa2PKCt', 'NgCaMCa2PKCp', 'NgCaMCa3', 'NgpCaMCa3', 'NgCaMCa3PKCt', 'NgCaMCa3PKCp', 'NgCaMCa4', 'NgpCaMCa4', 'NgCaMCa4PKCt', 'NgCaMCa4PKCp', 'NgpPP1', 'NgpPP2A', 'NgpPP2BCaMCa4', 'PP1', 'PP2A', 'CaM', 'CaMCa2', 'CaMCa3', 'CaMCa4', 'PP2B', 'PP2BCaM', 'PP2BCaMCa2', 'PP2BCaMCa3', 'PP2BCaMCa4', 'PKC', 'PKCCa', 'PKCt', 'PKCp' ]

interptimes = [10000*i for i in range(0,int((max(tvec))/10000))]
if interptimes[0] < 0:
  interptimes = interptimes[1:]
interpDATA = []
for j in range(0,len(species)):
  interpDATA.append(mytools.interpolate(tvec,vecs[j],interptimes))
tcDATA = array([interptimes]+interpDATA)
maxDATA = c_[tvec,array(vecs).T].max(axis=0)
scipy.io.savemat(filename, {'DATA': tcDATA, 'maxDATA': maxDATA, 'headers': headers})
