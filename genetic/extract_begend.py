import sys
from pylab import *

genelist_unord = sys.argv[1].split(',')
genelist = unique(genelist_unord)
begs = [100000000 for g in genelist]
ends = [0 for g in genelist]
chrs = [0 for g in genelist] 
chrs2 = [0 for g in genelist] 

input_file = open('ref_GRCh37.p13_top_level.gff3','r')
line = input_file.readline()
genelines = []
linechrs = []
linechrs2 = []
chrnow = 1
chr2now = 1
while line[0] == '#':
  line = input_file.readline()
while True:
  line = input_file.readline()
  if line.find('chromosome=') > -1:
    if line[line.find('chromosome=')+11] in ['0','1','2','3','4','5','6','7','8','9','X','Y']:
      chrnowstr = line[line.find('chromosome=')+11:line.find('chromosome=')+line[line.find('chromosome='):].find(';')]
      chrnow = int(chrnowstr) if chrnowstr != 'X' and chrnowstr != 'Y' else (23 if chrnowstr == 'X' else 24)
      #print('chrnow = '+str(chrnow)+', chrnowstr = '+chrnowstr)
  if line.find('chromosome ') > -1:
    if line[line.find('chromosome ')+11] in ['0','1','2','3','4','5','6','7','8','9','X','Y']:
      chr2nowstr = line[line.find('chromosome ')+11:line.find('chromosome ')+11+min(2,line[line.find('chromosome ')+11:].find(' '))]
      if chr2nowstr != '':
        chr2now = int(chr2nowstr) if chr2nowstr != 'X' and chr2nowstr != 'Y' else (23 if chr2nowstr == 'X' else 24)
        #print('chr2now = '+str(chr2now)+', chr2nowstr = '+chr2nowstr)


  if len(line) == 0:
    break
  splitted = line.split('\t')
  if len(splitted) < 9:
    continue
  entryType = splitted[2]
  if entryType == 'gene':
    genelines.append(line)
    linechrs.append(chrnow)
    linechrs2.append(chr2now)
    dataField = splitted[8]
    splitted2 = dataField.split(';')
    geneEntry = ''
    for isplitted in range(0,len(splitted2)):
      if splitted2[isplitted][0:5] == 'gene=':
        geneEntry = splitted2[isplitted][5:]

    igeneFound = -1
    if geneEntry in genelist:
      for igene in range(0,len(genelist)):
        if genelist[igene] == geneEntry:
          igeneFound = igene
          begs[igene] = int(splitted[3])
          ends[igene] = int(splitted[4])
          chrs[igene] = int(chrnow)
          chrs2[igene] = int(chr2now)
          break
      
input_file.close()
#print('chrnow = '+str(chrnow)+', chr2now = '+str(chr2now))
begs_syn = begs[:] #synonymous genes
ends_syn = ends[:] #synonymous genes
realnames = genelist[:]
chrs_syn = chrs[:]
chrs2_syn = chrs2[:]

for igene in range(0,len(genelist)):
  if ends[igene] == 0:
    for iline in range(0,len(genelines)):
      line = genelines[iline]
      splitted = line.split('\t')
      dataField = splitted[8]
      splitted2 = dataField.split(';')
      geneEntry = ''
      synonyms = ''
      for isplitted in range(0,len(splitted2)):
        if splitted2[isplitted][0:5] == 'gene=':
          geneEntry = splitted2[isplitted][5:]
        elif splitted2[isplitted][0:13] == 'gene_synonym=':
          synonyms = splitted2[isplitted][13:]
      splitted3 = synonyms.split(',')
      if genelist[igene] in splitted3:
        if ends_syn[igene] > 0:
          print('#Overwriting '+realnames[igene]+' for '+genelist[igene])
        realnames[igene] = geneEntry
        begs_syn[igene] = int(splitted[3])
        ends_syn[igene] = int(splitted[4])
        chrs_syn[igene] = linechrs[iline]
        chrs2_syn[igene] = linechrs2[iline]

#print('Only main gene names:')
#for i in range(0,len(genelist)):
#  print("'"+genelist[i]+"',"+str(chrs[i])+","+str(chrs2[i])+","+str(begs[i])+","+str(ends[i]))
#
#print('Using synonyms:')
#for i in range(0,len(genelist)):
#  print("'"+genelist[i]+"' ('"+realnames[i]+"'),"+str(chrs_syn[i])+","+str(chrs2_syn[i])+","+str(begs_syn[i])+","+str(ends_syn[i]))

print("DATA = [")
for i in range(0,len(genelist)):
  print("['"+realnames[i]+"',"+str(chrs_syn[i])+","+str(begs_syn[i])+","+str(ends_syn[i])+"],")
print("]")

#print(str(begs))
#print(str(ends))
