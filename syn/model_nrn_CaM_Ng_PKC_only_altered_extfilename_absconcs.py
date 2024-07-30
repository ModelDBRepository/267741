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

initvalues = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0]

species = ['Ca', 'Ng', 'NgCaM', 'Ngp', 'NgpCaM', 'NgPKCt', 'NgCaMPKCt', 'NgPKCp', 'NgCaMPKCp', 'NgCaMCa2', 'NgpCaMCa2', 'NgCaMCa2PKCt', 'NgCaMCa2PKCp', 'NgCaMCa3', 'NgpCaMCa3', 'NgCaMCa3PKCt', 'NgCaMCa3PKCp', 'NgCaMCa4', 'NgpCaMCa4', 'NgCaMCa4PKCt', 'NgCaMCa4PKCp', 'CaM', 'CaMCa2', 'CaMCa3', 'CaMCa4', 'PKC', 'PKCCa', 'PKCt', 'PKCp']

print("my_volume = "+str(my_volume)+" l ?= "+str(dend.L*(dend.diam/2)**2*3.14159265358)+" um3")
if len(initfile) > 3 and initfile != 'None':
  DATA_init = scipy.io.loadmat(initfile)
  for mykey in list(DATA_init.keys()):
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
ks = [1.0]*71
ks[0]   = 17000.0            # CaM + Ca*2 <-> CaMCa2 (forward)
ks[1]   = 0.035              # CaM + Ca*2 <-> CaMCa2 (backward)
ks[2]   = 14.0               # CaMCa2 + Ca <-> CaMCa3 (forward)
ks[3]   = 0.228              # CaMCa2 + Ca <-> CaMCa3 (backward)
ks[4]   = 26.0               # CaMCa3 + Ca <-> CaMCa4 (forward)
ks[5]   = 0.064              # CaMCa3 + Ca <-> CaMCa4 (backward)
ks[6]   = 15806.7            # NgCaM + Ca*2 <-> NgCaMCa2 (forward)
ks[7]   = 1.72117            # NgCaM + Ca*2 <-> NgCaMCa2 (backward)
ks[8]   = 14.0               # NgCaMCa2 + Ca <-> NgCaMCa3 (forward)
ks[9]   = 0.228              # NgCaMCa2 + Ca <-> NgCaMCa3 (backward)
ks[10]  = 26.0               # NgCaMCa3 + Ca <-> NgCaMCa4 (forward)
ks[11]  = 0.064              # NgCaMCa3 + Ca <-> NgCaMCa4 (backward)
ks[12]  = 15806.7            # NgpCaM + Ca*2 <-> NgpCaMCa2 (forward)
ks[13]  = 1.72117            # NgpCaM + Ca*2 <-> NgpCaMCa2 (backward)
ks[14]  = 14.0               # NgpCaMCa2 + Ca <-> NgpCaMCa3 (forward)
ks[15]  = 0.228              # NgpCaMCa2 + Ca <-> NgpCaMCa3 (backward)
ks[16]  = 26.0               # NgpCaMCa3 + Ca <-> NgpCaMCa4 (forward)
ks[17]  = 0.064              # NgpCaMCa3 + Ca <-> NgpCaMCa4 (backward)
ks[18]  = 28.0               # CaM + Ng <-> NgCaM (forward)
ks[19]  = 0.036              # CaM + Ng <-> NgCaM (backward)
ks[20]  = 0.0                # CaM + Ngp <-> NgpCaM (forward)
ks[21]  = 0.036              # CaM + Ngp <-> NgpCaM (backward)
ks[22]  = 2.0                # CaMCa2 + Ng <-> NgCaMCa2 (forward)
ks[23]  = 0.136              # CaMCa2 + Ng <-> NgCaMCa2 (backward)
ks[24]  = 0.0                # CaMCa2 + Ngp <-> NgpCaMCa2 (forward)
ks[25]  = 0.136              # CaMCa2 + Ngp <-> NgpCaMCa2 (backward)
ks[26]  = 2.0                # CaMCa3 + Ng <-> NgCaMCa3 (forward)
ks[27]  = 0.136              # CaMCa3 + Ng <-> NgCaMCa3 (backward)
ks[28]  = 0.0                # CaMCa3 + Ngp <-> NgpCaMCa3 (forward)
ks[29]  = 0.136              # CaMCa3 + Ngp <-> NgpCaMCa3 (backward)
ks[30]  = 2.0                # CaMCa4 + Ng <-> NgCaMCa4 (forward)
ks[31]  = 0.136              # CaMCa4 + Ng <-> NgCaMCa4 (backward)
ks[32]  = 0.0                # CaMCa4 + Ngp <-> NgpCaMCa4 (forward)
ks[33]  = 0.136              # CaMCa4 + Ngp <-> NgpCaMCa4 (backward)
ks[34]  = 13.299999999999999 # Ca + PKC <-> PKCCa (forward)
ks[35]  = 0.05               # Ca + PKC <-> PKCCa (backward)
ks[36]  = 0.0678             # Ng + PKCt <-> NgPKCt (forward)
ks[37]  = 0.002              # Ng + PKCt <-> NgPKCt (backward)
ks[38]  = 1.775              # NgPKCt --> Ngp + PKCt (forward)
ks[39]  = 0.0678             # NgCaM + PKCt <-> NgCaMPKCt (forward)
ks[40]  = 0.002              # NgCaM + PKCt <-> NgCaMPKCt (backward)
ks[41]  = 0.0678             # NgCaMCa2 + PKCt <-> NgCaMCa2PKCt (forward)
ks[42]  = 0.002              # NgCaMCa2 + PKCt <-> NgCaMCa2PKCt (backward)
ks[43]  = 0.0678             # NgCaMCa3 + PKCt <-> NgCaMCa3PKCt (forward)
ks[44]  = 0.002              # NgCaMCa3 + PKCt <-> NgCaMCa3PKCt (backward)
ks[45]  = 0.0678             # NgCaMCa4 + PKCt <-> NgCaMCa4PKCt (forward)
ks[46]  = 0.002              # NgCaMCa4 + PKCt <-> NgCaMCa4PKCt (backward)
ks[47]  = 1.775              # NgCaMPKCt --> NgpCaM + PKCt (forward)
ks[48]  = 1.775              # NgCaMCa2PKCt --> NgpCaMCa2 + PKCt (forward)
ks[49]  = 1.775              # NgCaMCa3PKCt --> NgpCaMCa3 + PKCt (forward)
ks[50]  = 1.775              # NgCaMCa4PKCt --> NgpCaMCa4 + PKCt (forward)
ks[51]  = 0.0678             # Ng + PKCp <-> NgPKCp (forward)
ks[52]  = 0.002              # Ng + PKCp <-> NgPKCp (backward)
ks[53]  = 1.775              # NgPKCp --> Ngp + PKCp (forward)
ks[54]  = 0.0678             # NgCaM + PKCp <-> NgCaMPKCp (forward)
ks[55]  = 0.002              # NgCaM + PKCp <-> NgCaMPKCp (backward)
ks[56]  = 0.0678             # NgCaMCa2 + PKCp <-> NgCaMCa2PKCp (forward)
ks[57]  = 0.002              # NgCaMCa2 + PKCp <-> NgCaMCa2PKCp (backward)
ks[58]  = 0.0678             # NgCaMCa3 + PKCp <-> NgCaMCa3PKCp (forward)
ks[59]  = 0.002              # NgCaMCa3 + PKCp <-> NgCaMCa3PKCp (backward)
ks[60]  = 0.0678             # NgCaMCa4 + PKCp <-> NgCaMCa4PKCp (forward)
ks[61]  = 0.002              # NgCaMCa4 + PKCp <-> NgCaMCa4PKCp (backward)
ks[62]  = 1.775              # NgCaMPKCp --> NgpCaM + PKCp (forward)
ks[63]  = 1.775              # NgCaMCa2PKCp --> NgpCaMCa2 + PKCp (forward)
ks[64]  = 1.775              # NgCaMCa3PKCp --> NgpCaMCa3 + PKCp (forward)
ks[65]  = 1.775              # NgCaMCa4PKCp --> NgpCaMCa4 + PKCp (forward)
ks[66]  = 2.5e-06            # Ngp --> Ng (forward)
ks[67]  = 2.5e-06            # NgpCaM --> NgCaM (forward)
ks[68]  = 2.5e-06            # NgpCaMCa2 --> NgCaMCa2 (forward)
ks[69]  = 2.5e-06            # NgpCaMCa3 --> NgCaMCa3 (forward)
ks[70]  = 2.5e-06            # NgpCaMCa4 --> NgCaMCa4 (forward)

for ialteredk in range(0,len(alteredks)):
  ks[alteredks[ialteredk]] = alteredk_factors[ialteredk]*ks[alteredks[ialteredk]]
reaction000 = rxd.Reaction(specs[21] + specs[0]*2 != specs[22], ks[0], ks[1])
reaction001 = rxd.Reaction(specs[22] + specs[0] != specs[23], ks[2], ks[3])
reaction002 = rxd.Reaction(specs[23] + specs[0] != specs[24], ks[4], ks[5])
reaction003 = rxd.Reaction(specs[2] + specs[0]*2 != specs[9], ks[6], ks[7])
reaction004 = rxd.Reaction(specs[9] + specs[0] != specs[13], ks[8], ks[9])
reaction005 = rxd.Reaction(specs[13] + specs[0] != specs[17], ks[10], ks[11])
reaction006 = rxd.Reaction(specs[4] + specs[0]*2 != specs[10], ks[12], ks[13])
reaction007 = rxd.Reaction(specs[10] + specs[0] != specs[14], ks[14], ks[15])
reaction008 = rxd.Reaction(specs[14] + specs[0] != specs[18], ks[16], ks[17])
reaction009 = rxd.Reaction(specs[21] + specs[1] != specs[2], ks[18], ks[19])
reaction010 = rxd.Reaction(specs[21] + specs[3] != specs[4], ks[20], ks[21])
reaction011 = rxd.Reaction(specs[22] + specs[1] != specs[9], ks[22], ks[23])
reaction012 = rxd.Reaction(specs[22] + specs[3] != specs[10], ks[24], ks[25])
reaction013 = rxd.Reaction(specs[23] + specs[1] != specs[13], ks[26], ks[27])
reaction014 = rxd.Reaction(specs[23] + specs[3] != specs[14], ks[28], ks[29])
reaction015 = rxd.Reaction(specs[24] + specs[1] != specs[17], ks[30], ks[31])
reaction016 = rxd.Reaction(specs[24] + specs[3] != specs[18], ks[32], ks[33])
reaction017 = rxd.Reaction(specs[0] + specs[25] != specs[26], ks[34], ks[35])
reaction018 = rxd.Reaction(specs[1] + specs[27] != specs[5], ks[36], ks[37])
reaction019 = rxd.Reaction(specs[5] > specs[3] + specs[27], ks[38])
reaction020 = rxd.Reaction(specs[2] + specs[27] != specs[6], ks[39], ks[40])
reaction021 = rxd.Reaction(specs[9] + specs[27] != specs[11], ks[41], ks[42])
reaction022 = rxd.Reaction(specs[13] + specs[27] != specs[15], ks[43], ks[44])
reaction023 = rxd.Reaction(specs[17] + specs[27] != specs[19], ks[45], ks[46])
reaction024 = rxd.Reaction(specs[6] > specs[4] + specs[27], ks[47])
reaction025 = rxd.Reaction(specs[11] > specs[10] + specs[27], ks[48])
reaction026 = rxd.Reaction(specs[15] > specs[14] + specs[27], ks[49])
reaction027 = rxd.Reaction(specs[19] > specs[18] + specs[27], ks[50])
reaction028 = rxd.Reaction(specs[1] + specs[28] != specs[7], ks[51], ks[52])
reaction029 = rxd.Reaction(specs[7] > specs[3] + specs[28], ks[53])
reaction030 = rxd.Reaction(specs[2] + specs[28] != specs[8], ks[54], ks[55])
reaction031 = rxd.Reaction(specs[9] + specs[28] != specs[12], ks[56], ks[57])
reaction032 = rxd.Reaction(specs[13] + specs[28] != specs[16], ks[58], ks[59])
reaction033 = rxd.Reaction(specs[17] + specs[28] != specs[20], ks[60], ks[61])
reaction034 = rxd.Reaction(specs[8] > specs[4] + specs[28], ks[62])
reaction035 = rxd.Reaction(specs[12] > specs[10] + specs[28], ks[63])
reaction036 = rxd.Reaction(specs[16] > specs[14] + specs[28], ks[64])
reaction037 = rxd.Reaction(specs[20] > specs[18] + specs[28], ks[65])
reaction038 = rxd.Reaction(specs[3] > specs[1], ks[66])
reaction039 = rxd.Reaction(specs[4] > specs[2], ks[67])
reaction040 = rxd.Reaction(specs[10] > specs[9], ks[68])
reaction041 = rxd.Reaction(specs[14] > specs[13], ks[69])
reaction042 = rxd.Reaction(specs[18] > specs[17], ks[70])
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

headers = [ 'tvec', 'Ca', 'Ng', 'NgCaM', 'Ngp', 'NgpCaM', 'NgPKCt', 'NgCaMPKCt', 'NgPKCp', 'NgCaMPKCp', 'NgCaMCa2', 'NgpCaMCa2', 'NgCaMCa2PKCt', 'NgCaMCa2PKCp', 'NgCaMCa3', 'NgpCaMCa3', 'NgCaMCa3PKCt', 'NgCaMCa3PKCp', 'NgCaMCa4', 'NgpCaMCa4', 'NgCaMCa4PKCt', 'NgCaMCa4PKCp', 'CaM', 'CaMCa2', 'CaMCa3', 'CaMCa4', 'PKC', 'PKCCa', 'PKCt', 'PKCp' ]

interptimes = [10000*i for i in range(0,int((max(tvec))/10000))]
if interptimes[0] < 0:
  interptimes = interptimes[1:]
interpDATA = []
for j in range(0,len(species)):
  interpDATA.append(mytools.interpolate(tvec,vecs[j],interptimes))
tcDATA = array([interptimes]+interpDATA)
maxDATA = c_[tvec,array(vecs).T].max(axis=0)
scipy.io.savemat(filename, {'DATA': tcDATA, 'maxDATA': maxDATA, 'headers': headers})
