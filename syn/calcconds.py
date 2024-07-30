import matplotlib
matplotlib.use('Agg')
from pylab import *
import scipy.io
import sys
import itertools
from os.path import exists
import mytools

def calcconds(filename, filename_nrn):
  #from sh getspecies_calcconds.sh|sed 's/, \]/\]/g'
  species = [ [ ['GluR1_S831', 'GluR1_S845_S831', 'GluR1_S831_PKAc', 'GluR1_S845_S831_PP1', 'GluR1_S831_PP1', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_S831_PP2B'] ],
              [ ['GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_S831_PP2B'] ],
              [ ['GluR1_S845', 'GluR1_S845_S831', 'GluR1_S845_CKCam', 'GluR1_S845_CKpCam', 'GluR1_S845_CKp', 'GluR1_S845_PKCtCa', 'GluR1_S845_PKCtAACa', 'GluR1_S845_PKCtDAGCa', 'GluR1_S845_PKCtAADAGCa', 'GluR1_S845_PP1', 'GluR1_S845_S831_PP1', 'GluR1_S845_PP2B', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S845', 'GluR1_memb_S845_S831', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_CKp', 'GluR1_memb_S845_PKCtCa', 'GluR1_memb_S845_PKCtAACa', 'GluR1_memb_S845_PKCtDAGCa', 'GluR1_memb_S845_PKCtAADAGCa', 'GluR1_memb_S845_PP1', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S845_PP2B', 'GluR1_memb_S845_S831_PP2B'] ],
              [ ['GluR1_memb', 'GluR1_memb_S845', 'GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_PKAc', 'GluR1_memb_CKCam', 'GluR1_memb_CKpCam', 'GluR1_memb_CKp', 'GluR1_memb_PKCtCa', 'GluR1_memb_PKCtAACa', 'GluR1_memb_PKCtDAGCa', 'GluR1_memb_PKCtAADAGCa', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_CKp', 'GluR1_memb_S845_PKCtCa', 'GluR1_memb_S845_PKCtAACa', 'GluR1_memb_S845_PKCtDAGCa', 'GluR1_memb_S845_PKCtAADAGCa', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_PP1', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_PP2B', 'GluR1_memb_S845_S831_PP2B'] ],
              [ ['GluR2_S880', 'GluR2_S880_PP2A', 'GluR2_memb_S880', 'GluR2_memb_S880_PP2A'] ],
              [ ['GluR2_memb', 'GluR2_memb_PKCtCa', 'GluR2_memb_PKCtAACa', 'GluR2_memb_PKCtDAGCa', 'GluR2_memb_PKCtAADAGCa', 'GluR2_memb_S880', 'GluR2_memb_S880_PP2A'] ] ]

  conds_hom1 = [12.4, 18.9]
  conds_hom2 = 2.2
  conds_het = 2.5
  Nskip = 1

  mesh_input_file = open('mesh_general.out','r')
  mesh_firstline = mesh_input_file.readline()
  mesh_secondline = mesh_input_file.readline()
  mesh_values = mesh_secondline.split()
  my_volume = float(mesh_values[-2])*1e-15 #litres
  mesh_input_file.close()

  MAT = {}
  assert exists(filename)
  MAT = scipy.io.loadmat(filename)

  DATA_all = MAT['DATA']
  header_strs = MAT['headers']
  for i in range(0,len(header_strs)):
    first_space = header_strs[i].find(' ')
    if first_space > -1:
      header_strs[i] = header_strs[i][:first_space]
  inddict = {}
  for iheader in range(4,len(header_strs)):
    inddict[header_strs[iheader]] = iheader-4

  DATANRN_all = {}
  assert exists(filename_nrn)
  DATANRN_all_all = scipy.io.loadmat(filename_nrn)
  for ikey in range(0,len(DATANRN_all_all['headers'])):
    if DATANRN_all_all['headers'][ikey].find(' ') > -1:
      mykey = DATANRN_all_all['headers'][ikey][0:DATANRN_all_all['headers'][ikey].find(' ')]
    else:
      mykey = DATANRN_all_all['headers'][ikey]
    DATANRN_all[mykey] = DATANRN_all_all['DATA'][ikey]

  if len(MAT) > 0:
    times = [500000/(DATA_all.shape[0]-1)*i for i in range(0,DATA_all.shape[0])]
  if len(DATANRN_all) > 0:
    times_nrn = DATANRN_all['tvec']

  TCs_all = []
  TCsN_all = []
  TCs_nrn_all = []
  TCsN_nrn_all = []
  for iax in range(0,len(species)):
    for ispecgroup in range(0,len(species[iax])):
      specgroup = species[iax][ispecgroup]
      if len(MAT) > 0:
        mytimecourse = zeros(DATA_all[:,0].shape[0])
      if len(DATANRN_all) > 0:
        mytimecourse_nrn = zeros(times_nrn.shape[0])
      if type(specgroup) is not list:
        specgroup = [specgroup]
      for ispec in range(0,len(specgroup)):
        specfactor = 1.0
        if len(MAT) > 0:
          mytimecourse = mytimecourse + specfactor*DATA_all[:,inddict[specgroup[ispec]]]
        if len(DATANRN_all) > 0:
          mytimecourse_nrn = mytimecourse_nrn + DATANRN_all[specgroup[ispec]]

      factor = 1.0/6.022e23/my_volume*1e9
      nrnfactor = 1.0
      TCs_nrn_all.append(mytimecourse_nrn[::Nskip]*1e6*nrnfactor)
      TCsN_nrn_all.append(mytimecourse_nrn[::Nskip]*1e6*nrnfactor/factor)
      TCs_all.append(mytimecourse[::Nskip]*factor)
      TCsN_all.append(mytimecourse[::Nskip])

  ENhom1_np_nrn = (TCsN_nrn_all[3] + TCsN_nrn_all[5])/4.0 * (TCsN_nrn_all[3]-TCsN_nrn_all[1])**4/(TCsN_nrn_all[3] + TCsN_nrn_all[5])**4                       
  ENhom1_p_nrn = (TCsN_nrn_all[3] + TCsN_nrn_all[5])/4.0 * (TCsN_nrn_all[3]**4 - (TCsN_nrn_all[3]-TCsN_nrn_all[1])**4)/(TCsN_nrn_all[3] + TCsN_nrn_all[5])**4 
  ENhom2_nrn = (TCsN_nrn_all[3] + TCsN_nrn_all[5])/4.0 * (TCsN_nrn_all[5]/(TCsN_nrn_all[3] + TCsN_nrn_all[5]))**4
  ENhet_nrn = (TCsN_nrn_all[3] + TCsN_nrn_all[5])/4.0 * (1 - (TCsN_nrn_all[3]/(TCsN_nrn_all[3] + TCsN_nrn_all[5]))**4 - (TCsN_nrn_all[5]/(TCsN_nrn_all[3] + TCsN_nrn_all[5]))**4)
  Egtot_nrn = ENhom1_np_nrn*conds_hom1[0] + ENhom1_p_nrn*conds_hom1[1] + ENhom2_nrn*conds_hom2 + ENhet_nrn*conds_het

  ENhom1_np = (TCsN_all[3] + TCsN_all[5])/4.0 * (TCsN_all[3]-TCsN_all[1])**4/(TCsN_all[3] + TCsN_all[5])**4                                #Number of complexes times the probability of a complex consisting of 4 non-phos GluR1s
  ENhom1_p = (TCsN_all[3] + TCsN_all[5])/4.0 * (TCsN_all[3]**4 - (TCsN_all[3]-TCsN_all[1])**4)/(TCsN_all[3] + TCsN_all[5])**4              #The same, but use prob. of having 4 GluR1s, minus the cases where all of them are non-phos
  ENhom2 = (TCsN_all[3] + TCsN_all[5])/4.0 * (TCsN_all[5]/(TCsN_all[3] + TCsN_all[5]))**4
  ENhet = (TCsN_all[3] + TCsN_all[5])/4.0 * (1 - (TCsN_all[3]/(TCsN_all[3] + TCsN_all[5]))**4 - (TCsN_all[5]/(TCsN_all[3] + TCsN_all[5]))**4)
  Egtot = ENhom1_np*conds_hom1[0] + ENhom1_p*conds_hom1[1] + ENhom2*conds_hom2 + ENhet*conds_het

  return [Egtot, Egtot_nrn, times, times_nrn]





def calcconds_nrn(filename_nrn):
  #from sh getspecies_calcconds.sh|sed 's/, \]/\]/g'
  species = [ [ ['GluR1_S831', 'GluR1_S845_S831', 'GluR1_S831_PKAc', 'GluR1_S845_S831_PP1', 'GluR1_S831_PP1', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_S831_PP2B'] ],
              [ ['GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_S831_PP2B'] ],
              [ ['GluR1_S845', 'GluR1_S845_S831', 'GluR1_S845_CKCam', 'GluR1_S845_CKpCam', 'GluR1_S845_CKp', 'GluR1_S845_PKCtCa', 'GluR1_S845_PKCtAACa', 'GluR1_S845_PKCtDAGCa', 'GluR1_S845_PKCtAADAGCa', 'GluR1_S845_PP1', 'GluR1_S845_S831_PP1', 'GluR1_S845_PP2B', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S845', 'GluR1_memb_S845_S831', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_CKp', 'GluR1_memb_S845_PKCtCa', 'GluR1_memb_S845_PKCtAACa', 'GluR1_memb_S845_PKCtDAGCa', 'GluR1_memb_S845_PKCtAADAGCa', 'GluR1_memb_S845_PP1', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S845_PP2B', 'GluR1_memb_S845_S831_PP2B'] ],
              [ ['GluR1_memb', 'GluR1_memb_S845', 'GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_PKAc', 'GluR1_memb_CKCam', 'GluR1_memb_CKpCam', 'GluR1_memb_CKp', 'GluR1_memb_PKCtCa', 'GluR1_memb_PKCtAACa', 'GluR1_memb_PKCtDAGCa', 'GluR1_memb_PKCtAADAGCa', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_CKp', 'GluR1_memb_S845_PKCtCa', 'GluR1_memb_S845_PKCtAACa', 'GluR1_memb_S845_PKCtDAGCa', 'GluR1_memb_S845_PKCtAADAGCa', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_PP1', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_PP2B', 'GluR1_memb_S845_S831_PP2B'] ],
              [ ['GluR2_S880', 'GluR2_S880_PP2A', 'GluR2_memb_S880', 'GluR2_memb_S880_PP2A'] ],
              [ ['GluR2_memb', 'GluR2_memb_PKCtCa', 'GluR2_memb_PKCtAACa', 'GluR2_memb_PKCtDAGCa', 'GluR2_memb_PKCtAADAGCa', 'GluR2_memb_S880', 'GluR2_memb_S880_PP2A'] ] ]

  conds_hom1 = [12.4, 18.9]
  conds_hom2 = 2.2
  conds_het = 2.5
  Nskip = 1

  mesh_input_file = open('mesh_general.out','r')
  mesh_firstline = mesh_input_file.readline()
  mesh_secondline = mesh_input_file.readline()
  mesh_values = mesh_secondline.split()
  my_volume = float(mesh_values[-2])*1e-15 #litres
  mesh_input_file.close()

  DATANRN_all = {}
  assert exists(filename_nrn)
  DATANRN_all_all = scipy.io.loadmat(filename_nrn)
  for ikey in range(0,len(DATANRN_all_all['headers'])):
    if DATANRN_all_all['headers'][ikey].find(' ') > -1:
      mykey = DATANRN_all_all['headers'][ikey][0:DATANRN_all_all['headers'][ikey].find(' ')]
    else:
      mykey = DATANRN_all_all['headers'][ikey]
    DATANRN_all[mykey] = DATANRN_all_all['DATA'][ikey]

  if len(DATANRN_all) > 0:
    times_nrn = DATANRN_all['tvec']

  TCs_nrn_all = []
  TCsN_nrn_all = []
  for iax in range(0,len(species)):
    for ispecgroup in range(0,len(species[iax])):
      specgroup = species[iax][ispecgroup]
      if len(DATANRN_all) > 0:
        mytimecourse_nrn = zeros(times_nrn.shape[0])
      if type(specgroup) is not list:
        specgroup = [specgroup]
      for ispec in range(0,len(specgroup)):
        specfactor = 1.0
        if len(DATANRN_all) > 0:
          mytimecourse_nrn = mytimecourse_nrn + DATANRN_all[specgroup[ispec]]

      factor = 1.0/6.022e23/my_volume*1e9
      nrnfactor = 1.0
      TCs_nrn_all.append(mytimecourse_nrn[::Nskip]*1e6*nrnfactor)
      TCsN_nrn_all.append(mytimecourse_nrn[::Nskip]*1e6*nrnfactor/factor)

  ENhom1_np_nrn = (TCsN_nrn_all[3] + TCsN_nrn_all[5])/4.0 * (TCsN_nrn_all[3]-TCsN_nrn_all[1])**4/(TCsN_nrn_all[3] + TCsN_nrn_all[5])**4                       
  ENhom1_p_nrn = (TCsN_nrn_all[3] + TCsN_nrn_all[5])/4.0 * (TCsN_nrn_all[3]**4 - (TCsN_nrn_all[3]-TCsN_nrn_all[1])**4)/(TCsN_nrn_all[3] + TCsN_nrn_all[5])**4 
  ENhom2_nrn = (TCsN_nrn_all[3] + TCsN_nrn_all[5])/4.0 * (TCsN_nrn_all[5]/(TCsN_nrn_all[3] + TCsN_nrn_all[5]))**4
  ENhet_nrn = (TCsN_nrn_all[3] + TCsN_nrn_all[5])/4.0 * (1 - (TCsN_nrn_all[3]/(TCsN_nrn_all[3] + TCsN_nrn_all[5]))**4 - (TCsN_nrn_all[5]/(TCsN_nrn_all[3] + TCsN_nrn_all[5]))**4)
  Egtot_nrn = ENhom1_np_nrn*conds_hom1[0] + ENhom1_p_nrn*conds_hom1[1] + ENhom2_nrn*conds_hom2 + ENhet_nrn*conds_het

  return Egtot_nrn, times_nrn



def calcconds_nrn_withcas(filename_nrn):
  #from sh getspecies_calcconds.sh|sed 's/, \]/\]/g'
  species = [ [ ['GluR1_S831', 'GluR1_S845_S831', 'GluR1_S831_PKAc', 'GluR1_S845_S831_PP1', 'GluR1_S831_PP1', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_S831_PP2B'] ],
              [ ['GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_S831_PP2B'] ],
              [ ['GluR1_S845', 'GluR1_S845_S831', 'GluR1_S845_CKCam', 'GluR1_S845_CKpCam', 'GluR1_S845_CKp', 'GluR1_S845_PKCtCa', 'GluR1_S845_PKCtAACa', 'GluR1_S845_PKCtDAGCa', 'GluR1_S845_PKCtAADAGCa', 'GluR1_S845_PP1', 'GluR1_S845_S831_PP1', 'GluR1_S845_PP2B', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S845', 'GluR1_memb_S845_S831', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_CKp', 'GluR1_memb_S845_PKCtCa', 'GluR1_memb_S845_PKCtAACa', 'GluR1_memb_S845_PKCtDAGCa', 'GluR1_memb_S845_PKCtAADAGCa', 'GluR1_memb_S845_PP1', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S845_PP2B', 'GluR1_memb_S845_S831_PP2B'] ],
              [ ['GluR1_memb', 'GluR1_memb_S845', 'GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_PKAc', 'GluR1_memb_CKCam', 'GluR1_memb_CKpCam', 'GluR1_memb_CKp', 'GluR1_memb_PKCtCa', 'GluR1_memb_PKCtAACa', 'GluR1_memb_PKCtDAGCa', 'GluR1_memb_PKCtAADAGCa', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_CKp', 'GluR1_memb_S845_PKCtCa', 'GluR1_memb_S845_PKCtAACa', 'GluR1_memb_S845_PKCtDAGCa', 'GluR1_memb_S845_PKCtAADAGCa', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_PP1', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_PP2B', 'GluR1_memb_S845_S831_PP2B'] ],
              [ ['GluR2_S880', 'GluR2_S880_PP2A', 'GluR2_memb_S880', 'GluR2_memb_S880_PP2A'] ],
              [ ['GluR2_memb', 'GluR2_memb_PKCtCa', 'GluR2_memb_PKCtAACa', 'GluR2_memb_PKCtDAGCa', 'GluR2_memb_PKCtAADAGCa', 'GluR2_memb_S880', 'GluR2_memb_S880_PP2A'] ] ]

  conds_hom1 = [12.4, 18.9]
  conds_hom2 = 2.2
  conds_het = 2.5
  Nskip = 1

  mesh_input_file = open('mesh_general.out','r')
  mesh_firstline = mesh_input_file.readline()
  mesh_secondline = mesh_input_file.readline()
  mesh_values = mesh_secondline.split()
  my_volume = float(mesh_values[-2])*1e-15 #litres
  mesh_input_file.close()

  DATANRN_all = {}
  assert exists(filename_nrn)
  DATANRN_all_all = scipy.io.loadmat(filename_nrn)
  for ikey in range(0,len(DATANRN_all_all['headers'])):
    if DATANRN_all_all['headers'][ikey].find(' ') > -1:
      mykey = DATANRN_all_all['headers'][ikey][0:DATANRN_all_all['headers'][ikey].find(' ')]
    else:
      mykey = DATANRN_all_all['headers'][ikey]
    DATANRN_all[mykey] = DATANRN_all_all['DATA'][ikey]

  if len(DATANRN_all) > 0:
    times_nrn = DATANRN_all['tvec']

  TCs_nrn_all = []
  TCsN_nrn_all = []
  for iax in range(0,len(species)):
    for ispecgroup in range(0,len(species[iax])):
      specgroup = species[iax][ispecgroup]
      if len(DATANRN_all) > 0:
        mytimecourse_nrn = zeros(times_nrn.shape[0])
      if type(specgroup) is not list:
        specgroup = [specgroup]
      for ispec in range(0,len(specgroup)):
        specfactor = 1.0
        #if len(specgroup[ispec]) > 24 and len(DATANRN_all) > 0:
        #  mytimecourse_nrn = mytimecourse_nrn + DATANRN_all[specgroup[ispec][:24]]
        #elif len(DATANRN_all) > 0:
        mytimecourse_nrn = mytimecourse_nrn + DATANRN_all[specgroup[ispec]]

      factor = 1.0/6.022e23/my_volume*1e9
      nrnfactor = 1.0
      TCs_nrn_all.append(mytimecourse_nrn[::Nskip]*1e6*nrnfactor)
      TCsN_nrn_all.append(mytimecourse_nrn[::Nskip]*1e6*nrnfactor/factor)

  ENhom1_np_nrn = (TCsN_nrn_all[3] + TCsN_nrn_all[5])/4.0 * (TCsN_nrn_all[3]-TCsN_nrn_all[1])**4/(TCsN_nrn_all[3] + TCsN_nrn_all[5])**4                       
  ENhom1_p_nrn = (TCsN_nrn_all[3] + TCsN_nrn_all[5])/4.0 * (TCsN_nrn_all[3]**4 - (TCsN_nrn_all[3]-TCsN_nrn_all[1])**4)/(TCsN_nrn_all[3] + TCsN_nrn_all[5])**4 
  ENhom2_nrn = (TCsN_nrn_all[3] + TCsN_nrn_all[5])/4.0 * (TCsN_nrn_all[5]/(TCsN_nrn_all[3] + TCsN_nrn_all[5]))**4
  ENhet_nrn = (TCsN_nrn_all[3] + TCsN_nrn_all[5])/4.0 * (1 - (TCsN_nrn_all[3]/(TCsN_nrn_all[3] + TCsN_nrn_all[5]))**4 - (TCsN_nrn_all[5]/(TCsN_nrn_all[3] + TCsN_nrn_all[5]))**4)
  Egtot_nrn = ENhom1_np_nrn*conds_hom1[0] + ENhom1_p_nrn*conds_hom1[1] + ENhom2_nrn*conds_hom2 + ENhet_nrn*conds_het

  return Egtot_nrn, times_nrn, DATANRN_all['Ca']


def calcconds_completePKC_nrn_withcas(filename_nrn):
  #from sh getspecies_calcconds.sh|sed 's/, \]/\]/g'
  species = [ [ ['GluR1_S831', 'GluR1_S845_S831', 'GluR1_S831_PKAc', 'GluR1_S845_S831_PP1', 'GluR1_S831_PP1', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_S831_PP2B'] ],
              [ ['GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_S831_PP2B'] ],
              [ ['GluR1_S845', 'GluR1_S845_S831', 'GluR1_S845_CKCam', 'GluR1_S845_CKpCam', 'GluR1_S845_CKp', 'GluR1_S845_PKCtCa', 'GluR1_S845_PKCtAACa', 'GluR1_S845_PKCtDAGCa', 'GluR1_S845_PKCtAADAGCa', 'GluR1_S845_PP1', 'GluR1_S845_S831_PP1', 'GluR1_S845_PP2B', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S845', 'GluR1_memb_S845_S831', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_CKp', 'GluR1_memb_S845_PKCtCa', 'GluR1_memb_S845_PKCtAACa', 'GluR1_memb_S845_PKCtDAGCa', 'GluR1_memb_S845_PKCtAADAGCa', 'GluR1_memb_S845_PP1', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S845_PP2B', 'GluR1_memb_S845_S831_PP2B'] ],
              [ ['GluR1_memb', 'GluR1_memb_S845', 'GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_PKAc', 'GluR1_memb_CKCam', 'GluR1_memb_CKpCam', 'GluR1_memb_CKp', 'GluR1_memb_PKCtCa', 'GluR1_memb_PKCtAACa', 'GluR1_memb_PKCtDAGCa', 'GluR1_memb_PKCtAADAGCa', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_CKp', 'GluR1_memb_S845_PKCtCa', 'GluR1_memb_S845_PKCtAACa', 'GluR1_memb_S845_PKCtDAGCa', 'GluR1_memb_S845_PKCtAADAGCa', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_PP1', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_PP2B', 'GluR1_memb_S845_S831_PP2B'] ],
              [ ['GluR2_S880', 'GluR2_S880_PP2A', 'GluR2_memb_S880', 'GluR2_memb_S880_PP2A'] ],
              [ ['GluR2_memb', 'GluR2_memb_PKCtCa', 'GluR2_memb_PKCtAACa', 'GluR2_memb_PKCtDAGCa', 'GluR2_memb_PKCtAADAGCa', 'GluR2_memb_S880', 'GluR2_memb_S880_PP2A'] ] ]

  conds_hom1 = [12.4, 18.9]
  conds_hom2 = 2.2
  conds_het = 2.5
  Nskip = 1

  mesh_input_file = open('mesh_general.out','r')
  mesh_firstline = mesh_input_file.readline()
  mesh_secondline = mesh_input_file.readline()
  mesh_values = mesh_secondline.split()
  my_volume = float(mesh_values[-2])*1e-15 #litres
  mesh_input_file.close()

  DATANRN_all = {}
  assert exists(filename_nrn)
  DATANRN_all_all = scipy.io.loadmat(filename_nrn)
  for ikey in range(0,len(DATANRN_all_all['headers'])):
    if DATANRN_all_all['headers'][ikey].find(' ') == -1:
      mykey = DATANRN_all_all['headers'][ikey]
    else:
      mykey = DATANRN_all_all['headers'][ikey][0:DATANRN_all_all['headers'][ikey].find(' ')]
    DATANRN_all[mykey] = DATANRN_all_all['DATA'][ikey]

  if len(DATANRN_all) > 0:
    times_nrn = DATANRN_all['tvec']

  TCs_nrn_all = []
  TCsN_nrn_all = []
  for iax in range(0,len(species)):
    for ispecgroup in range(0,len(species[iax])):
      specgroup = species[iax][ispecgroup]
      if len(DATANRN_all) > 0:
        mytimecourse_nrn = zeros(times_nrn.shape[0])
      if type(specgroup) is not list:
        specgroup = [specgroup]
      for ispec in range(0,len(specgroup)):
        specfactor = 1.0
        #if len(specgroup[ispec]) > 24 and len(DATANRN_all) > 0:
        #  mytimecourse_nrn = mytimecourse_nrn + DATANRN_all[specgroup[ispec][:24]]
        #elif len(DATANRN_all) > 0:
        mytimecourse_nrn = mytimecourse_nrn + DATANRN_all[specgroup[ispec]]

      factor = 1.0/6.022e23/my_volume*1e9
      nrnfactor = 1.0
      TCs_nrn_all.append(mytimecourse_nrn[::Nskip]*1e6*nrnfactor)
      TCsN_nrn_all.append(mytimecourse_nrn[::Nskip]*1e6*nrnfactor/factor)

  ENhom1_np_nrn = (TCsN_nrn_all[3] + TCsN_nrn_all[5])/4.0 * (TCsN_nrn_all[3]-TCsN_nrn_all[1])**4/(TCsN_nrn_all[3] + TCsN_nrn_all[5])**4                       
  ENhom1_p_nrn = (TCsN_nrn_all[3] + TCsN_nrn_all[5])/4.0 * (TCsN_nrn_all[3]**4 - (TCsN_nrn_all[3]-TCsN_nrn_all[1])**4)/(TCsN_nrn_all[3] + TCsN_nrn_all[5])**4 
  ENhom2_nrn = (TCsN_nrn_all[3] + TCsN_nrn_all[5])/4.0 * (TCsN_nrn_all[5]/(TCsN_nrn_all[3] + TCsN_nrn_all[5]))**4
  ENhet_nrn = (TCsN_nrn_all[3] + TCsN_nrn_all[5])/4.0 * (1 - (TCsN_nrn_all[3]/(TCsN_nrn_all[3] + TCsN_nrn_all[5]))**4 - (TCsN_nrn_all[5]/(TCsN_nrn_all[3] + TCsN_nrn_all[5]))**4)
  Egtot_nrn = ENhom1_np_nrn*conds_hom1[0] + ENhom1_p_nrn*conds_hom1[1] + ENhom2_nrn*conds_hom2 + ENhet_nrn*conds_het

  return Egtot_nrn, times_nrn, DATANRN_all['Ca']

