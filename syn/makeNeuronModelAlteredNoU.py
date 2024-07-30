from os.path import exists

bodyname = 'nrn_altered_noU'
output_file = open('model_'+str(bodyname)+'.py','w')

### Write the header including the initiatlization and reading of arguments etc.
output_file.write("""from neuron import h, rxd
from pylab import *
from matplotlib import pyplot
import scipy.io
import time
import re
import mytools

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
tolerance = 1e-6
Ca_input_onset = 800
Ca_input_N     = 100
Ca_input_freq  = 100
Ca_input_dur   = 0.005
Ca_input_flux  = 600.0
L_input_flux   = 2.0
Glu_input_flux = 50.0
ACh_input_flux = 2.0
Ntrains        = 1
trainT = 3000
initfile = ''
addition = ''
blocked = []
blockeds = []
block_factor = 1.0
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
  Ca_input_onset = float(sys.argv[3])
if len(sys.argv) > 4:
  Ca_input_N     = int(sys.argv[4])
if len(sys.argv) > 5:
  Ca_input_freq  = float(sys.argv[5])
if len(sys.argv) > 6:
  Ca_input_dur   = float(sys.argv[6])
if len(sys.argv) > 7:
  Ca_input_flux  = float(sys.argv[7])
if len(sys.argv) > 8:
  L_input_flux   = float(sys.argv[8])
if len(sys.argv) > 9:
  Glu_input_flux = float(sys.argv[9])
if len(sys.argv) > 10:
  ACh_input_flux = float(sys.argv[10])
if len(sys.argv) > 11:
  Ntrains  = int(float(sys.argv[11]))
if len(sys.argv) > 12:
  trainT  = float(sys.argv[12])
if len(sys.argv) > 13:
  initfile = sys.argv[13]
if len(sys.argv) > 14:
  blocked = sys.argv[14]
  blockeds = blocked.split(',')
if len(sys.argv) > 15:
  block_factor = sys.argv[15]
  block_factors = [float(x) for x in block_factor.split(',')]
if type(blocked) is not list:
  addition = '_'+blocked+'x'+str(block_factor)
if len(sys.argv) > 16:
  alteredk = sys.argv[16]
  alteredks = [int(x) for x in alteredk.split(',')]
if len(sys.argv) > 17:
  alteredk_factor = sys.argv[17]
  alteredk_factors = [float(x) for x in alteredk_factor.split(',')]
if type(alteredk) is not list:
  addition = addition+'_k'+alteredk+'x'+str(alteredk_factor)
if len(sys.argv) > 18:
  toltochange = 'Gs,Calbin,CalbinC,PMCA,PMCACa,PP2B,Gi,GiaGDP,Gibg,Gsbg,GsaGDP,CaOut,CaOutLeak,Leak,PDE1,Ng,NgCaM,NgpCaM,CaM,CaMCa2,CaMCa4,fixedbuffer,fixedbufferCa,DGL,CaDGL,'+sys.argv[18]
  tolstochange = toltochange.split(',')
if len(sys.argv) > 19:
  tolchange = '3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,2,2,1,1,'+sys.argv[19]
  tolschange = [float(x) for x in tolchange.split(',')]
""")

### Read the initial concentrations
IC_file = open('IC_singlecompartment.xml')
firstline = IC_file.readline()
line = IC_file.readline()
species_all = []
initvals_all = []
while len(line) > 0:
  if line.find('NanoMolarity') > -1 and (line.find('<!--') < 0 or line.find('NanoMolarity') < line.find('<!--')):
    species = line[line.find('"')+1:line.find('"')+1+line[line.find('"')+1:].find('"')]
    if species[0] > '0' and species[0] <= '9':
      species = '_'+species
    restline = line[line.find('"')+1+line[line.find('"')+1:].find('"')+1:]
    value = float(restline[restline.find('"')+1:restline.find('"')+1+restline[restline.find('"')+1:].find('"')])
    print('spec'+str(len(species_all))+' (spec_'+species+') = rxd.Species(cyt, name=\''+species+'\', initial='+str(value*1e-6))
    initvals_all.append(value*1e-6)
    species_all.append(species)
  line = IC_file.readline()
IC_file.close()

### Write the initial concentrations to 'initvalues' and the corresponding species to 'species'
output_file.write('\ninitvalues = [')
for ispec in range(0,len(initvals_all)-1):
  output_file.write(str(initvals_all[ispec])+', ')
output_file.write(str(initvals_all[-1])+']\n')
output_file.write('\nspecies = [')
for ispec in range(0,len(species_all)-1):
  if species_all[ispec][0] == '_':
    output_file.write('\''+species_all[ispec][1:]+'\', ')
  else:
    output_file.write('\''+species_all[ispec]+'\', ')
if species_all[-1][0] == '_':
  output_file.write('\''+species_all[-1][1:]+'\']\n')
else:
  output_file.write('\''+species_all[-1]+'\']\n')

### Change initvalues if needed
output_file.write("""
print("my_volume = "+str(my_volume)+" l ?= "+str(dend.L*(dend.diam/2)**2*3.14159265358)+" um3")
if len(initfile) > 3 and initfile != 'None':
  DATA_init = scipy.io.loadmat(initfile)
  for ispec in range(0,len(species)):
    initvalues[ispec] = DATA_init['DATA'][1+ispec,-1]
    if species[ispec] not in DATA_init['headers'][1+ispec]:
      print("Warning: mismatch in DATA_init, ispec = "+str(ispec))
tolscales = [1.0 for i in range(0,len(species))]
for ispec in range(0,len(species)):
  for iblock in range(0,len(blockeds)):
    if re.match(blockeds[iblock]+'$',species[ispec]):
      initvalues[ispec] = block_factors[iblock]*initvalues[ispec]
  for itol in range(0,len(tolstochange)):
    if tolstochange[itol] == species[ispec]:
      tolscales[ispec] = tolscales[ispec]*10**(-tolschange[itol])
specs = []
for ispec in range(0,len(species)):
  specs.append(rxd.Species(cyt, name='spec'+str(ispec), charge=0, initial=initvalues[ispec], atolscale=tolscales[ispec]))
  if tolscales[ispec] != 1.0:
    print('spec'+str(ispec)+' ('+species[ispec]+') atolscale = '+str(tolscales[ispec]))
""")

output_file.write('Ca_flux_rate = rxd.Parameter(cyt, initial=0)\n')
output_file.write('L_flux_rate = rxd.Parameter(cyt, initial=0)\n')
output_file.write('Glu_flux_rate = rxd.Parameter(cyt, initial=0)\n')
output_file.write('ACh_flux_rate = rxd.Parameter(cyt, initial=0)\n')


Reacs = []
Reactants = []
Reactants_int = []
Reactant_powers = []
Reactant_ns = []
Products = []
Products_int = []
Product_powers = []
Product_ns = []
Forwardrates = []
Backwardrates = []
reac_file = open('Reactions.xml')
firstline = reac_file.readline()
line = reac_file.readline()

ireaction = 0
ireaction_str = '000'
reactionRatesAll = []       # Includes all reaction rates, forward and backward
reactionForwardsAll = []    # Includes string 'forward' for rates that are forward and 'backward' for rates that are backward
ireactionsAll = []          # Includes indices for reactions. There will be more reaction rates than reactions (as some reactions contain forw. and backw. rates), this is the mapping between them
reactionLinesAll = []       # Includes the lines to be output to the output_file. These are saved first and written later as the reaction rates are allowed to be changed later
reactionLinesAllShort = []  # Includes short versions of the lines in reactionLinesAll. These will be in comments - with real specie names; it helps the programmer to know which rate affects which reaction

### Read the reactions
while len(line) > 0:
  if line.find('Reaction name') > -1 and line[0:4] != '<!--':
    newline = reac_file.readline()

    Reactants_this = []
    Products_this = []
    Reactants_this_int = []
    Products_this_int = []
    Forwardrate = 0
    Backwardrate = 0
    Powers = []
    PowersProd = []
    ns = []
    nsProd = []
    while newline.find('</Reaction>') == -1:
      if newline.find('Reactant specieID') > -1 or newline.find('Reactant  specieID') > -1:
        Reactants_this.append(newline[newline.find('"')+1:newline.find('"')+newline[newline.find('"')+1:].find('"')+1])
        if Reactants_this[-1][0] > '0' and Reactants_this[-1][0] <= '9':
          Reactants_this[-1] = '_'+Reactants_this[-1]
        if newline.find('power') > -1:
          newline2 = newline[newline.find('power'):]
          Powers.append(int(newline2[newline2.find('"')+1:newline2.find('"')+newline2[newline2.find('"')+1:].find('"')+1]))
        else:
          Powers.append(1)
        if newline.find('n=') > -1:
          newline2 = newline[newline.find('n='):]
          ns.append(int(newline2[newline2.find('"')+1:newline2.find('"')+newline2[newline2.find('"')+1:].find('"')+1]))
        else:
          ns.append(Powers[-1])

      if newline.find('Product specieID') > -1 or newline.find('Product  specieID') > -1:
        Products_this.append(newline[newline.find('"')+1:newline.find('"')+newline[newline.find('"')+1:].find('"')+1])
        if Products_this[-1][0] > '0' and Products_this[-1][0] <= '9':
          Products_this[-1] = '_'+Products_this[-1]
        if newline.find('power') > -1:
          newline2 = newline[newline.find('power'):]
          PowersProd.append(int(newline2[newline2.find('"')+1:newline2.find('"')+newline2[newline2.find('"')+1:].find('"')+1]))
        else:
          PowersProd.append(1)
        if newline.find('n=') > -1:
          newline2 = newline[newline.find('n='):]
          nsProd.append(int(newline2[newline2.find('"')+1:newline2.find('"')+newline2[newline2.find('"')+1:].find('"')+1]))
        else:
          nsProd.append(PowersProd[-1])
      if newline.find('forwardRate') > -1:
        newline2 = newline[newline.find('>')+1:]
        Forwardrate = float(newline2[:newline2.find('<')])
      if newline.find('reverseRate') > -1:
        newline2 = newline[newline.find('>')+1:]
        Backwardrate = float(newline2[:newline2.find('<')])
      newline = reac_file.readline()
    #print("while ended: newline = "+newline)
    if len(Products_this) == 0:
      Products_this.append(Reactants_this[0])
      PowersProd.append(0)
      nsProd.append(0)

    for i in range(0,len(Reactants_this)):
      foundone = 0
      for j in range(0,len(species_all)):
        if species_all[j] == Reactants_this[i]:
          foundone = 1
          break
      if not foundone:
        for j in range(0,len(species_all)):
          if species_all[j] == Reactants_this[i][1:]:
            foundone = 1
            print("found _")
            break
      if not foundone:
        print("Error: "+Reactants_this[i]+" not found in species!!!!")
        time.sleep(10)
      Reactants_this_int.append(j)

    for i in range(0,len(Products_this)):
      foundone = 0
      for j in range(0,len(species_all)):
        if species_all[j] == Products_this[i]:
          foundone = 1
          break
      if not foundone:
        for j in range(0,len(species_all)):
          if species_all[j] == Products_this[i][1:]:
            foundone = 1
            print("found _")
            break
      if not foundone:
        print("Error: "+Products_this[i]+" not found in species!!!!")
        time.sleep(10)
      Products_this_int.append(j)

    Forwardrates.append(Forwardrate)
    Backwardrates.append(Backwardrate)
    Reactants.append(Reactants_this[:])
    Reactants_int.append(Reactants_this_int[:])
    Products.append(Products_this[:])
    Products_int.append(Products_this_int[:])
    Reactant_powers.append(Powers[:])
    Product_powers.append(PowersProd[:])
    Reactant_ns.append(ns[:])
    Product_ns.append(nsProd[:])
    #print("Reaction #"+str(len(Reactants))+" added")

    #Go through the reactants and form the left side of the reaction arrow. Take into account that if n=something, this has to be included
    Reac_txt = ''
    Reac_txt_short = ''
    for i in range(0,len(Reactants_this)):
      Reac_txt = Reac_txt + 'specs['+str(Reactants_this_int[i])+']'
      Reac_txt_short = Reac_txt_short + Reactants_this[i]
      if ns[i] != 1:
        Reac_txt = Reac_txt + "*" + str(ns[i])
        Reac_txt_short = Reac_txt_short + "*" + str(ns[i])
      if i < len(Reactants_this) - 1:
        Reac_txt = Reac_txt + " + "
        Reac_txt_short = Reac_txt_short + " + "
    #Reac_txt = Reac_txt + " = "

    #Check if all Powers are the same as all ns. If not, the reaction rate has to be made custom, based on Powers instead of ns
    if ns == Powers and nsProd == PowersProd:
      #rate_txt = ', '+str(Forwardrate*(1e6)**(sum([int(x) for x in Powers])-1))
      rate_txt = ', ks['+str(len(reactionRatesAll))+']'
      ireactionsAll.append(ireaction)
      reactionRatesAll.append(Forwardrate*(1e6)**(sum([int(x) for x in Powers])-1))
      reactionForwardsAll.append('forward')
    else:
      spec_prod_txt = ''
      for i in range(0,len(Reactants_this)):
        spec_prod_txt = spec_prod_txt+'specs['+str(Reactants_this_int[i])+']'
        if int(Powers[i]) != 1:
          spec_prod_txt = spec_prod_txt+'**'+str(Powers[i])
        if i < len(Reactants_this)-1:
          spec_prod_txt = spec_prod_txt+'*'
      #rate_txt = ', '+str(Forwardrate*(1e6)**(sum([int(x) for x in Powers])-1))+'*'+spec_prod_txt
      rate_txt = ', ks['+str(len(reactionRatesAll))+']*'+spec_prod_txt
      ireactionsAll.append(ireaction)
      reactionRatesAll.append(Forwardrate*(1e6)**(sum([int(x) for x in Powers])-1))
      reactionForwardsAll.append('forward')
          
    if Backwardrate == 0:
      Reac_txt = Reac_txt + ", "
      Reac_txt_short = Reac_txt_short + " --> "
    else:
      Reac_txt = Reac_txt + ", " #Change this to != if you are using Python 3
      Reac_txt_short = Reac_txt_short + " <-> "
      #Check if all Powers are the same as all ns. If not, the reaction rate has to be made custom, based on PowersProd instead of nsProd
      if ns == Powers and nsProd == PowersProd:
        rate_txt = rate_txt+', ks['+str(len(reactionRatesAll))+']'
        ireactionsAll.append(ireaction)
        reactionRatesAll.append(Backwardrate*(1e6)**(sum([int(x) for x in PowersProd])-1))
        reactionForwardsAll.append('backward')
      else:
        spec_prod_txt = ''
        for i in range(0,len(Products_this)):
          spec_prod_txt = spec_prod_txt+'specs['+str(Products_this_int[i])+']'
          if int(PowersProd[i]) != 1:
            spec_prod_txt = spec_prod_txt+'**'+str(PowersProd[i])
          if i < len(Products_this)-1:
            spec_prod_txt = spec_prod_txt+'*'
        rate_txt = rate_txt+', ks['+str(len(reactionRatesAll))+']*'+spec_prod_txt
        ireactionsAll.append(ireaction)
        reactionRatesAll.append(Backwardrate*(1e6)**(sum([int(x) for x in PowersProd])-1))
        reactionForwardsAll.append('backward')

    #Go through the products and form the right side of the reaction arrow. Take into account that if n=something, this has to be included
    for i in range(0,len(Products_this)):
      Reac_txt = Reac_txt + 'specs['+str(Products_this_int[i])+']'
      Reac_txt_short = Reac_txt_short + Products_this[i]
      if nsProd[i] != 1:
        Reac_txt = Reac_txt + "*" + str(nsProd[i])
        Reac_txt_short = Reac_txt_short + "*" + str(nsProd[i])
      if i < len(Products_this) - 1:
        Reac_txt = Reac_txt + " + "
        Reac_txt_short = Reac_txt_short + " + "

    #If ns different from Powers, remember to assign custom_dynamics=True (otherwise the mass-action reaction rate(s) will be multiplied by the given rate term(s))
    if ns == Powers and nsProd == PowersProd:
      print('reaction'+ireaction_str+' = rxd.Reaction('+Reac_txt+rate_txt+')')
      #output_file.write('reaction'+ireaction_str+' = rxd.Reaction('+Reac_txt+rate_txt+')\n')
      reactionLinesAll.append('reaction'+ireaction_str+' = rxd.Reaction('+Reac_txt+rate_txt+')\n')
    else:
      print('reaction'+ireaction_str+' = rxd.Reaction('+Reac_txt+rate_txt+', custom_dynamics=True)')
      #output_file.write('reaction'+ireaction_str+' = rxd.Reaction('+Reac_txt+rate_txt+', custom_dynamics=True)\n')
      reactionLinesAll.append('reaction'+ireaction_str+' = rxd.Reaction('+Reac_txt+rate_txt+', custom_dynamics=True)\n')
    reactionLinesAllShort.append(Reac_txt_short)
    ireaction = ireaction + 1
    ireaction_str = str(ireaction)
    if ireaction < 10:
      ireaction_str = '0'+ireaction_str
    if ireaction < 100:
      ireaction_str = '0'+ireaction_str
    Reacs.append(Reac_txt)
  line = reac_file.readline()
    
reac_file.close()

### Write the reaction rates
output_file.write('ks = [1.0]*'+str(len(reactionRatesAll))+'\n')
for ireacRate in range(0,len(reactionRatesAll)):
  output_file.write('ks['+str(ireacRate)+']'+' '*(3-len(str(ireacRate)))+' = '+str(reactionRatesAll[ireacRate])+' '*(max([len(str(x)) for x in reactionRatesAll])-len(str(reactionRatesAll[ireacRate])))+
                    ' # '+reactionLinesAllShort[ireactionsAll[ireacRate]]+' ('+reactionForwardsAll[ireacRate]+')\n')

### Change reaction rates if needed
output_file.write("""
for ialteredk in range(0,len(alteredks)):
  ks[alteredks[ialteredk]] = alteredk_factors[ialteredk]*ks[alteredks[ialteredk]]
""")

### Write the reactions
for iline in range(0,len(reactionLinesAll)):
  output_file.write(reactionLinesAll[iline])

ispec_Ca = -1
ispec_L = -1
ispec_Glu = -1
ispec_ACh = -1
for ispec in range(0,len(species_all)):
  if species_all[ispec] == 'Ca':
    ispec_Ca = ispec
  if species_all[ispec] == 'L':
    ispec_L = ispec
  if species_all[ispec] == 'Glu':
    ispec_Glu = ispec
  if species_all[ispec] == 'ACh':
    ispec_ACh = ispec

### Introduce the variables that determine the input to the spine
output_file.write("""
reaction_Ca_flux = rxd.Rate(specs["""+str(ispec_Ca)+"""], Ca_flux_rate) # Ca
reaction_L_flux = rxd.Rate(specs["""+str(ispec_L)+"""], L_flux_rate) # L
reaction_Glu_flux = rxd.Rate(specs["""+str(ispec_Glu)+"""], Glu_flux_rate) # Glu
reaction_ACh_flux = rxd.Rate(specs["""+str(ispec_ACh)+"""], ACh_flux_rate) # ACh
""")

### Record the time courses
output_file.write('vec_t = h.Vector()\n')
output_file.write("""
vecs = []
vec_t = h.Vector()
vec_t.record(h._ref_t)
for ispec in range(0,len(species)):
  vecs.append(h.Vector())
  vecs[ispec].record(specs[ispec].nodes(dend)(0.5)[0]._ref_concentration)
""")

output_file.write("""
cvode = h.CVode()
cvode.active(1)
hmax = cvode.maxstep(1000)
hmin = cvode.minstep(1e-10)
cvode.atol(tolerance)

h.finitialize(-65)
def set_param(param, val):
    param.nodes.value = val
    h.cvode.re_init()

### Set on and off the inputs to the spine
T = 1000./Ca_input_freq
tnow = 0
for itrain in range(0,Ntrains):
    for istim in range(0,Ca_input_N):
      tnew = Ca_input_onset + istim*T + trainT*itrain
      h.cvode.event(tnew, lambda: set_param(Ca_flux_rate, Ca_input_flux/6.022e23/my_volume*1e3))
      h.cvode.event(tnew+Ca_input_dur, lambda: set_param(Ca_flux_rate, 0))
      h.cvode.event(tnew, lambda: set_param(L_flux_rate, L_input_flux/6.022e23/my_volume*1e3))
      h.cvode.event(tnew+Ca_input_dur, lambda: set_param(L_flux_rate, 0))
      h.cvode.event(tnew, lambda: set_param(Glu_flux_rate, Glu_input_flux/6.022e23/my_volume*1e3))
      h.cvode.event(tnew+Ca_input_dur, lambda: set_param(Glu_flux_rate, 0))
      h.cvode.event(tnew, lambda: set_param(ACh_flux_rate, ACh_input_flux/6.022e23/my_volume*1e3))
      h.cvode.event(tnew+Ca_input_dur, lambda: set_param(ACh_flux_rate, 0))
      tnow = tnew
timenow = time.time()
h.continuerun(Duration)
print("Simulation done in "+str(time.time()-timenow)+" seconds")
def isFlux(t):
  for itrain in range(0,Ntrains):
    for istim in range(0,Ca_input_N):
      tnew = Ca_input_onset + istim*T + trainT*itrain
      if t >= tnew and t < tnew+Ca_input_dur:
        return 1
  return 0
tvec = array(vec_t)
minDT_nonFlux = 20.0
minDT_Flux = 1.0
lastt = -inf
itvec2 = []
for it in range(0,len(tvec)):
  if tvec[it] - lastt > minDT_nonFlux or (isFlux(tvec[it]) and tvec[it] - lastt > minDT_Flux):
    itvec2.append(it)
    lastt = tvec[it]
""")

output_file.write("\nheaders = [ ")
output_file.write("'tvec', ")
for ispec in range(0,len(species_all)-1):
  if species_all[ispec][0] == '_':
    output_file.write("'"+species_all[ispec][1:]+"', ")
  else:
    output_file.write("'"+species_all[ispec]+"', ")
if species_all[-1][0] == '_':
  output_file.write("'"+species_all[-1][1:]+"' ]\n")
else:
  output_file.write("'"+species_all[-1]+"' ]\n")
output_file.write("""
myonset = Ca_input_onset
if myonset > max(tvec):
  myonset = 0
interptimes = [myonset + 10000*i for i in range(-1,int((max(tvec)-myonset)/10000))]
if interptimes[0] < 0:
  interptimes = interptimes[1:]
interpDATA = []
for j in range(0,len(species)):
  interpDATA.append(mytools.interpolate(tvec,vecs[j],interptimes))
tcDATA = array([interptimes]+interpDATA)
maxDATA = c_[tvec,array(vecs).T].max(axis=0)
filename = 'nrn_tstop'+str(Duration)+'_tol'+str(tolerance)+addition+'_onset'+str(Ca_input_onset)+'_n'+str(Ca_input_N)+'_freq'+str(Ca_input_freq)+'_dur'+str(Ca_input_dur)+'_flux'+str(Ca_input_flux)+'_Lflux'+str(L_input_flux)+'_Gluflux'+str(Glu_input_flux)+'_AChflux'+str(ACh_input_flux)+'_Ntrains'+str(Ntrains)+'_trainT'+str(trainT)+'.mat'
toBeRemovedIfNecessary = ['_tol1e-06','_tstop15000000','3560000_600000','_Ninputs1','_pulseamp5.0']
for i in range(0,len(toBeRemovedIfNecessary)):
  if len(filename) > 254:
    filename = filename.replace(toBeRemovedIfNecessary[i],'')
scipy.io.savemat(filename, {'DATA': tcDATA, 'maxDATA': maxDATA, 'headers': headers})
""")
output_file.close()


