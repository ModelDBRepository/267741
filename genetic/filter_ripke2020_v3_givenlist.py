import sys
import csv
import scipy.io

genelist = 'genes_new'
if len(sys.argv) > 1:
  genelist = sys.argv[1]

print("Running filtered_"+genelist+".py")
exec(open("filtered_"+genelist+".py").read())

filename = 'PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv'
input_file = open(filename,'r')
line = input_file.readline()
geneNamesAll = [DATA[i][0] for i in range(0,len(DATA))]
chrsAll = [DATA[i][1] for i in range(0,len(DATA))]
begsAll = [DATA[i][2] for i in range(0,len(DATA))]
endsAll = [DATA[i][3] for i in range(0,len(DATA))]
iline = 0
igenesFound = []
ilinesFound = []
rssFound = []
chrsFound = []
pvalsFound = []
while line[0] == '#' and len(line) > 0:
  line = input_file.readline()
while len(line) > 0:
  if iline % 10000 == 9999:
    print(str(iline+1)+"/7500000")
  line = input_file.readline()
  fields = line.split('\t')
  if len(fields) < 3:
    continue
  if len(fields[0]) == 0 or len(fields[1]) == 0 or len(fields[2]) == 0:
    print('iline = '+str(iline)+'; line = '+line+' has problems')
    continue
  for igene in range(0,len(geneNamesAll)):
    if chrsAll[igene] == int(fields[0]) and begsAll[igene] <= int(fields[2]) <= endsAll[igene]:
       igenesFound.append(igene)
       ilinesFound.append(iline)
       rssFound.append(fields[1])
       chrsFound.append(fields[0])
       pvalsFound.append(fields[10])
  iline = iline + 1
input_file.close()
print("igenesFound = "+str(igenesFound))
for iigene in range(0,len(igenesFound)):
  igene = igenesFound[iigene]
  if float(pvalsFound[iigene]) < 1e-6:
    print("  Gene="+str(geneNamesAll[igene])+" line="+str(ilinesFound[iigene])+" rs="+str(rssFound[iigene])+" chr="+str(chrsFound[iigene])+" pval="+str(pvalsFound[iigene]))
scipy.io.savemat('filtered_ripke2020_v3_'+genelist+'.mat',{'igenes': igenesFound, 'ilines': ilinesFound, 'rss': rssFound, 'chrs': chrsFound, 'pvals': pvalsFound, 'geneNamesAll': geneNamesAll}) #The results should be more or less the same
