# cp extract_commonmind_ACC_data_sklearn_allgenes_params_samples.py extract_commonmind_data_sklearn_allgenes_params_samples.py
from pylab import *
from sklearn.linear_model import LogisticRegression, LinearRegression
import statsmodels.api as sm
import pandas as pd
from scipy import stats

brainArea = 'DLPFC' #AnCg, nAcc
normalise = 1
if len(sys.argv) > 1:
  brainArea = sys.argv[1]
if len(sys.argv) > 2:
  normalise = int(sys.argv[2])
  
interestingGeneNames = ['CACNA1A', 'CACNA1B', 'CACNA1C', 'CACNA1D', 'CACNA1E', 'CACNA1F', 'CACNA1G', 'CACNA1H', 'CACNA1I', 'CACNA1S', 'CACNB1', 'CACNB2', 'CACNB3', 'CACNB4', 'CACNA2D1', 'CACNA2D2', 'CACNA2D3', 'CACNA2D4', 'CACNG1', 'CACNG2', 'CACNG3', 'CACNG4', 'CACNG5', 'CACNG6', 'CACNG7', 'CACNG8', 'SCN1A', 'SCN2A', 'SCN3A', 'SCN4A', 'SCN5A', 'SCN7A', 'SCN8A', 'SCN9A', 'SCN10A', 'SCN11A', 'SCN1B', 'SCN2B', 'SCN3B', 'SCN4B', 'KCNN1', 'KCNN2', 'KCNN3', 'KCNN4', 'ATP2A1', 'ATP2A2', 'ATP2A3', 'ATP2B1', 'ATP2B2', 'ATP2B3', 'ATP2B4', 'ATP2C1', 'ATP2C2', 'NRGN', 'KCNA1', 'KCNA2', 'KCNA3', 'KCNA4', 'KCNA5', 'KCNA6', 'KCNA7', 'KCNA10', 'KCNB1', 'KCNB2', 'KCNC1', 'KCNC2', 'KCNC3', 'KCNC4', 'KCND1', 'KCND2', 'KCND3', 'KCNQ1', 'KCNQ2', 'KCNQ3', 'KCNQ4', 'KCNQ5', 'KCNH1', 'KCNH2', 'KCNH3', 'KCNH4', 'KCNH5', 'KCNH6', 'KCNH7', 'KCNH8', 'KCNF1', 'KCNG1', 'KCNG2', 'KCNG3', 'KCNG4', 'KCNJ1', 'KCNJ2', 'KCNJ3', 'KCNJ4', 'KCNJ5', 'KCNJ6', 'KCNJ8', 'KCNJ9', 'KCNJ10', 'KCNJ11', 'KCNJ12', 'KCNJ13', 'KCNJ14', 'KCNJ15', 'KCNJ16', 'KCNJ17', 'KCNK1', 'KCNK2', 'KCNK3', 'KCNK4', 'KCNK5', 'KCNK6', 'KCNK7', 'KCNK9', 'KCNK10', 'KCNK12', 'KCNK13', 'KCNK15', 'KCNK16', 'KCNK17', 'KCNK18', 'KCNS1', 'KCNS2', 'KCNS3', 'KCNT1', 'KCNT2', 'KCNU1', 'KCNMA1', 'KCNMB1', 'KCNMB2', 'KCNMB3', 'KCNMB4', 'KCNV1', 'KCNV2', 'KCNAB1', 'KCNAB2', 'KCNAB3', 'KCNE1', 'KCNE2', 'KCNE3', 'KCNE4', 'KCNE1L', 'KCNIP1', 'KCNIP2', 'KCNIP3', 'KCNIP4', 'HCN1', 'HCN2', 'HCN3', 'HCN4', 'NKAIN1', 'NKAIN2', 'NKAIN3', 'NKAIN4', 'NKAIN1P1', 'NKAIN1P2', 'ITPR1', 'ITPR2', 'ITPR3', 'ITPRIP', 'SLC8A1', 'SLC8A2', 'SLC8A3', 'SLC8B1', 'SLC24A1', 'SLC24A2', 'SLC24A3', 'SLC24A4', 'SLC24A5', 'SLC32A1', 'SLC44A1', 'SLC44A2', 'SLC44A3', 'SLC44A4', 'SLC44A5', 'SLC44A3-AS1', 'SLC44A3P1', 'SLC12A1', 'SLC12A2', 'SLC12A3', 'SLC12A4', 'SLC12A5', 'SLC12A6', 'SLC12A7', 'SLC12A8', 'SLC12A9', 'SLC12A5-AS1', 'SLC12A9-AS1', 'CATSPER1', 'CATSPER2', 'CATSPER2P1', 'CATSPER3', 'CATSPER4', 'CATSPERB', 'CATSPERD', 'CATSPERE', 'CATSPERG', 'CATSPERZ', 'TRPC1', 'TRPC2', 'TRPC3', 'TRPC4', 'TRPC5', 'TRPC6', 'TRPC7', 'TRPC8', 'TRPC4AP', 'TRPC4OS', 'TRPC7-AS1', 'TRPC7-AS2', 'TRPM1', 'TRPM2', 'TRPM3', 'TRPM4', 'TRPM5', 'TRPM6', 'TRPM7', 'TRPM8', 'TRPV1', 'TRPV2', 'TRPV3', 'TRPV4', 'TRPV5', 'TRPV6', 'GRIA1', 'GRIA2', 'GRIA3', 'GRIA4', 'GRID1', 'GRID2', 'GRID1-AS1', 'GRID2IP', 'GRIN1', 'GRIN2A', 'GRIN2B', 'GRIN2C', 'GRIN2D', 'GRIN3A', 'GRIN3B', 'GCOM2', 'GABRA1', 'GABRA2', 'GABRA3', 'GABRA4', 'GABRA5', 'GABRA6', 'GABRB1', 'GABRB2', 'GABRB3', 'GABRD', 'GABRE', 'GABRG1', 'GABRG2', 'GABRG3', 'GABRG3-AS1', 'GABRP', 'GABRQ', 'GABRR1', 'GABRR2', 'GABRR3', 'GABBR1', 'GABBR2', 'CHRNA1', 'CHRNA2', 'CHRNA3', 'CHRNA4', 'CHRNA5', 'CHRNA6', 'CHRNA7', 'CHRNA9', 'CHRNA10', 'CHRNB1', 'CHRNB2', 'CHRNB3', 'CHRNB4', 'CHRND', 'CHRNE', 'CHRNG', 'GRM1', 'GRM2', 'GRM3', 'GRM4', 'GRM5', 'PRKACA', 'PRKACB', 'PRKACG', 'PRKAR1A', 'PRKAR1B', 'PRKAR2A', 'PRKAR2B', 'AKAP1', 'AKAP2', 'AKAP3', 'AKAP4', 'AKAP5', 'AKAP6', 'AKAP7', 'AKAP8', 'AKAP9', 'AKAP10', 'AKAP11', 'AKAP12', 'AKAP13', 'AKAP14', 'AKAP15', 'AKAP17A', 'AKAP17ABP', 'SPHKAP', 'PRKCA', 'PRKCB', 'PRKCD', 'PRKCE', 'PRKCG', 'PRKCH', 'PRKCI', 'PRKCQ', 'PRKCZ', 'PRKD1', 'PRKD2', 'PRKD3', 'PRKG1', 'PRKG2', 'CAMK2A', 'CAMK2B', 'CAMK2D', 'CAMK2G', 'ADCY1', 'ADCY2', 'ADCY3', 'ADCY4', 'ADCY5', 'ADCY6', 'ADCY7', 'ADCY8', 'ADCY9', 'ADCY10', 'PLCB1', 'PLCB2', 'PLCB3', 'PLCB4', 'PLCD1', 'PLCD3', 'PLCD4', 'PLCE1', 'PLCG1', 'PLCG2', 'PLCH1', 'PLCH2', 'PLCL1', 'PLCL2', 'CALM1', 'CALM2', 'CALM3', 'CALN1', 'GNAS', 'GNAS-AS1', 'GNAL', 'GNAQ', 'GNA11', 'GNA14', 'GNA14-AS1', 'GNA15', 'GNA12', 'GNA13', 'GNAI1', 'GNAI2', 'GNAI3', 'GNAO1', 'GNAZ', 'GNB1', 'GNB2', 'GNB3', 'GNB4', 'GNB5', 'GNB1L', 'GNG2', 'GNG3', 'GNG4', 'GNG5', 'GNG7', 'GNG8', 'GNG10', 'GNG11', 'GNG12', 'GNG13', 'GNG14', 'RGS1', 'RGS2', 'RGS3', 'RGS4', 'RGS5', 'RGS6', 'RGS7', 'RGS7BP', 'RGS8', 'RGS9', 'RGS9BP', 'RGS10', 'RGS11', 'RGS12', 'RGS13', 'RGS14', 'RGS16', 'RGS17', 'RGS18', 'RGS19', 'RGS20', 'RGS21', 'RGS22', 'GPSM1', 'GPSM2', 'GPSM3', 'GPS1', 'GPS2', 'GPS2P1', 'GPS2P2', 'CREB1', 'ATF4', 'CREB3', 'CREB5', 'CREB3L1', 'CREB3L2', 'CREB3L3', 'CREB3L4', 'CREBL2', 'ANP32A', 'ANP32B', 'ANP32C', 'ANP32D', 'ANP32E', 'ANP32A-IT1', 'ANP32AP1', 'ANP32BP2', 'ANP32BP3', 'ADRB1', 'ADRB2', 'ADRB3', 'ADRA1A', 'ADRA1B', 'ADRA1D', 'ADRA2A', 'ADRA2B', 'ADRA2C', 'DRD1', 'DRD2', 'DRD3', 'DRD4', 'DRD5', 'HTR1A', 'HTR1B', 'HTR1D', 'HTR1DP1', 'HTR1E', 'HTR1F', 'HTR2A', 'HTR2B', 'HTR2C', 'HTR3A', 'HTR3B', 'HTR3C', 'HTR3D', 'HTR3E', 'HTR4', 'HTR5A', 'HTR6', 'HTR7', 'OPRD1', 'OPRK1', 'OPRL1', 'OPRM1', 'OPRPN', 'CHRM1', 'CHRM2', 'CHRM3', 'CHRM4', 'CHRM5', 'CALB1', 'CALB2', 'S100G', 'RAPGEF3', 'ATP2B1', 'ATP2B2', 'ATP2B3', 'ATP2B4', 'SLC8A1', 'SLC8A2', 'SLC8A3', 'NRGN', 'PPP1R1A', 'PPP1R1B', 'PPP1R1C', 'PPP1CA', 'PPP1CB', 'PPP1CC', 'PPP1R2A', 'PPP1R2B', 'PPP1R2C', 'PPP1R3A', 'PPP1R3B', 'PPP1R3C', 'PPP1R3D', 'PPP1R3E', 'PPP1R3F', 'PPP1R3G', 'PPP1R7', 'PPP1R8', 'PPP1R9A', 'PPP1R9B', 'PPP1R10', 'PPP1R11', 'PPP1R12A', 'PPP1R12B', 'PPP1R12C', 'PPP1R13B', 'PPP1R13L', 'PPP1R14A', 'PPP1R14B', 'PPP1R14C', 'PPP1R14D', 'PPP1R15A', 'PPP1R15B', 'PPP1R16A', 'PPP1R16B', 'PHACTR1', 'PHACTR2', 'PHACTR3', 'PHACTR4', 'PDE1A', 'PDE1B', 'PDE1C', 'PDE2A', 'PDE3A', 'PDE3B', 'PDE4A', 'PDE4B', 'PDE4B-AS1', 'PDE4C', 'PDE4D', 'PDE4DIP', 'PDE4DIPP1', 'PDE5A', 'PDE6A', 'PDE6B', 'PDE6C', 'PDE6D', 'PDE6G', 'PDE6H', 'PDE7A', 'PDE7B', 'PDE8A', 'PDE8B', 'PDE9A', 'PDE10A', 'PDE11A', 'PPP2CA', 'PPP2CB', 'PPP2R1A', 'PPP2R1B', 'PPP2R2A', 'PPP2R2B', 'PPP2R2C', 'PPP2R2D', 'PPP2R3A', 'PPP2R3B', 'PPP2R3C', 'PPP2R4', 'PPP2R5A', 'PPP2R5B', 'PPP2R5C', 'PPP2R5D', 'PPP2R5E', 'PPP3CA', 'PPP3CB', 'PPP3CC', 'PPP3R1', 'PPP3R2', 'DGKA', 'DGKB', 'DGKD', 'DGKE', 'DGKG', 'DGKH', 'DGKI', 'DGKK', 'DGKQ', 'DGKZ', 'PLA2R1', 'PLA2G2A', 'PLA2G2C', 'PLA2G2D', 'PLA2G2E', 'PLA2G2F', 'PLA2G3', 'PLA2G4A', 'PLA2G4B', 'PLA2G4C', 'PLA2G4D', 'PLA2G4E', 'PLA2G4F', 'PLA2G5', 'PLA2G6', 'PLA2G7', 'PLA2G10', 'PLA2G12A', 'PLA2G12B', 'PLA2G15', 'PITPNM1', 'PITPNM2', 'PITPNM3', 'PITPNA', 'PITPNB', 'PITPNC1', 'NBEA', 'NBEAL1', 'NBEAL2', 'NBEAP1', 'NBEAP2', 'NBEAP3', 'NBEAP4', 'NBEAP5', 'NBEAP6', 'CSNK1A1', 'CSNK1A1L', 'CSNK1A1P1', 'CSNK1A1P2', 'CSNK1A1P3', 'CSNK1D', 'CSNK1E', 'CSNK1G1', 'CSNK1G2', 'CSNK1G3', 'CSNK1G2-AS1', 'CSNK1G2P1', 'CSNK2A1', 'CSNK2A2', 'CSNK2A3', 'CSNK2B', 'CSNKAIP', 'CDK5', 'CDK5R1', 'CDK5R2', 'CDK5RAP1', 'CDK5RAP2', 'CDK5RAP3', 'CDK5P1', 'CLK1', 'CLK2', 'CLK3', 'CLK4', 'CLK2P1', 'CLK3P1', 'CLK3P2', 'SCGN', 'NUCB1', 'NUCB2', 'NUCB1-AS1', 'AKT1', 'AKT2', 'AKT3', 'AKTIP', 'AKT1S1', 'CLTA', 'CLTB', 'CLTC', 'CLTCL1', 'CLTRN', 'SYNGAP1', 'DLG1', 'DLG2', 'DLG3', 'DLG4', 'DLG5', 'MPP', 'DLGAP5', 'DLGAP1', 'DLGAP2', 'DLGAP3', 'DLGAP4', 'DLGAP5', 'DLG1-AS1', 'DLG2-AS2', 'DLG3-AS1', 'DLG5-AS1', 'SRR', 'OPCML', 'OPCML-IT1', 'OPCML-IT2', 'NOS1', 'NOS2', 'NOS3', 'NOSIP', 'NOSTRIN', 'NOS1AP', 'NOS2P1', 'NOS2P2', 'NOS2P3', 'NOS2P4']
#interestingGeneNames = ['PPP3R1','CAMK2G','DGKI','GNAI2','RGS6','GNG7','GRM1','GRM5','SLC8A1','NRGN','PDE4B','PDE4D','PDE4DIP','AKAP6','PRKAR2A','PRKCB','PLCE1','PLCL2','PPP1R1A','PPP1R12C','PPP1R13B','PPP1R14C','PPP1R16B','PPP2R2B','PPP2R3A'] #the order should not affect the outcome
#TODO: make the right data correspond to right gene!

input_file = open('CommonMind_SCZ_metaData.txt','r')
line = input_file.readline()
samples = []
while len(line) > 0:
  line = input_file.readline()
  if len(line) < 10:
    #print("Short line")
    continue
  splitted = line.split('\t')
  samples.append([splitted[i] for i in [1,3,7,10,16]])
input_file.close()


#Remove XXY
for isample in range(len(samples)-1,-1,-1):
  if samples[isample][2] not in ['"XX"','"XY"']:
    samples.pop(isample)
  samples[isample][0] = samples[isample][0].replace('"','')

#input_file = open('CommonMind_SCZ_normalizedExpression.txt','r')
input_file = open('CIBERSORTxHiRes_PFC_Neurons_Window20.txt','r')
line = input_file.readline()
samples_expdata = line.split()[1:]
DATA = []
DATAigenes = []
lines_interesting = []

while len(line) > 0:
  line = input_file.readline()
  if len(line) < 2:
    continue
  splitted = line.split()
  if splitted[0] not in interestingGeneNames:
    #print('Omit '+splitted[0])
    continue
  if splitted[1] == 'NA':
    continue
  iinterestingGene = [i for i in range(0,len(interestingGeneNames)) if interestingGeneNames[i] == splitted[0]][0]
  lines_interesting.append(line)
  DATAigenes.append(iinterestingGene)
  DATA.append([float(x) for x in splitted[1:]])
input_file.close()

#Find in which order the gene data is in DATA and save the data in correct order into exprData
igeneord = []
genes_included = []
for i in range(0,len(interestingGeneNames)):
  for j in range(0,len(DATAigenes)):
    if DATAigenes[j] == i:
      igeneord.append(j)
      genes_included.append(interestingGeneNames[i])
      break

exprDataUnord = array(DATA)
exprData = exprDataUnord[igeneord,:]

genes = ['NRGN', 'SLC8A1', 'SLC8A3', 'GRM1', 'PRKCB', 'GNAI2', 'CALB2', 'PPP1CC', 'PDE4B', 'DGKI', 'DGKZ', 'CACNA1C', 'CACNA1D', 'CACNA1I', 'KCNB1', 'KCNQ3', 'HCN1', 'PRKAR2A', 'DGKH', 'GNAQ', 'CAMK2B', 'SCN9A', 'SCN1B', 'KCNK4', 'HCN1', 'DGKG', 'PPP1CA', 'PLA2G4A', 'PRKAR2A', 'PRKCA', 'SCN9A', 'CACNA1D', 'CACNA1C', 'KCND3']

betas = []
betas_sm = []
pvals = []
FDRs = []
FDRs_sm = []
isubjs_HC = [[i for i in range(0,len(samples_expdata)) if samples_expdata[i] == samples[isample][0]][0] for isample in range(0,len(samples)) if samples[isample][1] == '"Control"']
isubjs_SCZ = [[i for i in range(0,len(samples_expdata)) if samples_expdata[i] == samples[isample][0]][0] for isample in range(0,len(samples)) if samples[isample][1] == '"SCZ"']
print("isubjs_HC = "+str(isubjs_HC))
print("isubjs_SCZ = "+str(isubjs_SCZ))
meanExps = []
meanExps_HC = []
meanExps_SCZ = []
foundGeneNames = []
iexprs = []
igene_exprData = -1
for igene in range(0,len(interestingGeneNames)):
  if igene%10 == 0:
    print("igene = "+str(igene)+"/"+str(len(interestingGeneNames)))
  if igene not in DATAigenes:
    continue
  igene_exprData = igene_exprData + 1
  if interestingGeneNames[igene] not in genes: #This has to be done after "igene_exprData = igene_exprData + 1", otherwise correct numbering is lost
    continue
  if interestingGeneNames[igene] != genes_included[igene_exprData]:
    print('Error, check this')
    time.sleep(10)
  if interestingGeneNames[igene] != genes_included[igene_exprData]:
    print('Error, check this')
    time.sleep(10)
  X = []
  Y = []
  for isample in range(0,len(samples)):
    icols = [i for i in range(0,len(samples_expdata)) if samples_expdata[i] == samples[isample][0]]
    if len(icols) != 1:
      qwe
    isubj = icols[0]
    X.append([0 if samples[isample][1] == '"Control"' else (1 if samples[isample][1] == '"SCZ"' else 2),0 if samples[isample][2] == '"XX"' else 1, float(samples[isample][3]), float(samples[isample][4])])
    Y.append(exprData[igene_exprData,isubj])
  Xdf = pd.DataFrame(array(X),columns = ['Diagnosis','Sex','PMI','Age'])
  Y_HC_mean = mean([exprData[igene_exprData,i] for i in isubjs_HC])
  Ydf = pd.DataFrame(array(Y)/Y_HC_mean,columns = ['Expression'])

  if len(unique(Ydf['Expression'])) < 2:
    FDRs.append(nan)
    betas.append(array([[nan]*4]))
    betas_sm.append(array([[nan]*4]))
    pvals.append([nan,nan,nan,nan])
    meanExps.append(mean(Ydf['Expression']))
    meanExps_HC.append(mean(Ydf['Expression']))
    meanExps_SCZ.append(mean(Ydf['Expression']))
    continue

  #Do by statsmodels.api
  Xdf_sm = sm.add_constant(Xdf)
  log_reg_sm = sm.OLS(Ydf, Xdf_sm, missing='Drop')
  log_reg_sm.exog = sm.add_constant(log_reg_sm.exog, prepend=False)
  log_reg_sm_fit = log_reg_sm.fit()

  #Do by sklearn
  clf = LinearRegression().fit(Xdf, array(Ydf))
  predicted_category = clf.predict(Xdf)
  predicted_SCZ_correct = [predicted_category[i] and samples[i][1] == '"SCZ"' for i in range(0,len(samples))]
  predicted_SCZ_incorrect = [not predicted_category[i] and samples[i][1] == '"SCZ"' for i in range(0,len(samples))]
  FDRs.append(sum(predicted_SCZ_incorrect)/(sum(predicted_SCZ_incorrect) + sum(predicted_SCZ_correct)))
  if  abs(clf.intercept_ - log_reg_sm_fit.params[0])[0] > 0.01 or sum(abs(array(log_reg_sm_fit.params[1:])-clf.coef_[0])) > 0.01:
    print("clf.intercept = "+str(clf.intercept_)+", log_reg_sm_fit.params[0] = "+str(log_reg_sm_fit.params[0])+", array(log_reg_sm_fit.params[1:]) = "+str(array(log_reg_sm_fit.params[1:]))+", clf.coef_[0] = "+str(clf.coef_[0]))
  betas.append(clf.coef_)
  betas_sm.append(array(log_reg_sm_fit.params[1:]))
  pvals.append(log_reg_sm_fit.pvalues[:])
  #log_reg = LogisticRegression(random_state=0)
  #log_reg_fit = log_reg.fit(Xfd,Ydf)
  #print(log_reg_fit.summary())
  #betas.append(log_reg_fit.params[:])
  #pvals.append(log_reg_fit.pvalues[:])
  predicted_category = [z>0.5 for z in log_reg_sm.predict(log_reg_sm_fit.params[:],Xdf_sm)]
  #predicted_correct = [predicted_category[i] and samples[i][1] == '"SCZ"' or not predicted_category[i] and samples[i][1] == '"Control"' for i in range(0,len(samples))]
  #predicted_incorrect = [predicted_category[i] and samples[i][1] == '"Control"' or not predicted_category[i] and samples[i][1] == '"SCZ"' for i in range(0,len(samples))]
  predicted_SCZ_correct = [predicted_category[i] and samples[i][1] == '"SCZ"' for i in range(0,len(samples))]
  predicted_SCZ_incorrect = [not predicted_category[i] and samples[i][1] == '"SCZ"' for i in range(0,len(samples))]
  #predicted_HC_correct = [not predicted_category[i] and samples[i][1] == '"Control"' for i in range(0,len(samples))]
  #predicted_HC_incorrect = [predicted_category[i] and samples[i][1] == '"Control"' for i in range(0,len(samples))]
  FDRs_sm.append(sum(predicted_SCZ_incorrect)/(sum(predicted_SCZ_incorrect) + sum(predicted_SCZ_correct)))
  meanExps.append(mean(Ydf['Expression']))
  meanExps_HC.append(mean([exprData[igene_exprData,i] for i in isubjs_HC]))
  meanExps_SCZ.append(mean([exprData[igene_exprData,i] for i in isubjs_SCZ]))
  foundGeneNames.append(interestingGeneNames[igene])
  print(str(DATAigenes[igene_exprData]))
  iexprs.append(igene_exprData)


variables_genes = ['Ng', 'NCX', 'NCX', 'MGluR', 'PKC', 'Gi', 'Calbin*', 'PP1', 'PDE4', 'DAGK', 'DAGK', 'gCa_HVAbar_Ca_HVA', 'gCa_HVAbar_Ca_HVA', 'gCa_LVAstbar_Ca_LVAst', 'gK_Pstbar_K_Pst', 'gImbar_Im', 'gIhbar_Ih', 'PKA', 'DAGK', 'Gqabg', 'CK', 'gNaT*_tbar_NaT*_t', 'gNaT*_tbar_NaT*_t', 'g_pas', 'gIhbar_Ih', 'DAGK', 'PP1', 'PLA2', 'PKA', 'PKC', 'gNaT*_tbar_NaT*_t', 'gCa_HVAbar_Ca_HVA', 'gCa_HVAbar_Ca_HVA', 'gK_Pstbar_K_Pst']
#variables_genes = ['Ng', 'NCX', 'NCX', 'MGluR', 'PKC', 'Gi', 'Calbin*', 'PP1', 'PDE4', 'DAGK', 'DAGK', 'gCa_HVAbar_Ca_HVA', 'gCa_HVAbar_Ca_HVA', 'gCa_LVAstbar_Ca_LVAst', 'gK_Pstbar_K_Pst', 'gImbar_Im', 'gIhbar_Ih', 'PKA', 'DAGK', 'Gqabg', 'CK', 'gNap_Et2bar_Nap_Et2', 'gNaT*_tbar_NaT*_t', 'g_pas', 'gIhbar_Ih', 'DAGK', 'PP1', 'PLA2', 'PKA', 'PKC', 'gNap_Et2bar_Nap_Et2', 'gCa_HVAbar_Ca_HVA', 'gCa_HVAbar_Ca_HVA', 'gK_Pstbar_K_Pst']
variables = unique(variables_genes)
gene_to_variable_dict = {}
gene_to_igene_dict = {}
ivars_genes = []
for igene in range(0,len(genes)):
  gene_to_variable_dict[genes[igene]] = variables_genes[igene]
  gene_to_igene_dict[genes[igene]] = igene
  for ivar in range(0,len(variables)):
    if variables_genes[igene] == variables[ivar]:
      ivars_genes.append(ivar)

filename = 'genes_all_review_summary.csv'
input_file = open(filename,'r')
line = input_file.readline()
lines = []
PFCDATA = []
ACCDATA = []
genes2 = []
PFCpvals = []
ACCpvals = []
while len(line) > 0:
  line = input_file.readline()
  if len(line) < 2:
    continue
  splitted = line.split(';')
  lines.append(line)
  genes2.append(splitted[0])
  PFCDATA.append([float(x) for x in splitted[6:9]])
  ACCDATA.append([float(x) for x in splitted[9:12]])
  PFCpvals.append([float(x) for x in [splitted[12],splitted[14]]])
  ACCpvals.append([float(x) for x in [splitted[13],splitted[15]]])

input_file.close()
#ACCDATA
#list(zip(genes,ACCDATA))
#mean(exprData_all[0][:,isubjs_SCZ],axis=1)/mean(exprData_all[0][:,isubjs_HC],axis=1)


if len(genes) != len(genes2) or sum([genes2[i]!=genes[i] for i in range(0,len(genes))]) > 0:
    print('Check the gene indexing!')
    exit('Check the gene indexing!')
iGWASgenes = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
iPFConlygenes = [17,18,19,20,21,22,23,24]
iACConlygenes = [25,26,27,28,29,30,31,32,33]

iPFCgenes = iGWASgenes+iPFConlygenes
iACCgenes = iGWASgenes+iACConlygenes

for i in range(len(iPFCgenes)-1,-1,-1):
  if genes[iPFCgenes[i]] in [genes[iPFCgenes[j]] for j in range(0,i)]:
    iPFCgenes.pop(i)
for i in range(len(iACCgenes)-1,-1,-1):
  if genes[iACCgenes[i]] in [genes[iACCgenes[j]] for j in range(0,i)]:
    iACCgenes.pop(i)    

isImputed = 1
area = 'PFC'
exprCoeffs_all = []
exprData_all = []
igenes_all = []
isSignificants_all = []
isImputeds_all = []
geneSources_all = []
areas_all = []
for isSignificant in [0,1]:
  for geneSource in ['GWAS','CommonMind','Both']:
    myDATA = PFCDATA if area == 'PFC' else ACCDATA
    mypvals = PFCpvals if area == 'PFC' else ACCpvals
        
    iCommonMindgenes = iPFConlygenes if area == 'PFC' else iACConlygenes
    iBothgenes = iPFCgenes if area == 'PFC' else iACCgenes
    myigenes = iGWASgenes if geneSource == 'GWAS' else (iCommonMindgenes if geneSource == 'CommonMind' else iBothgenes)
        
    igenes = [i for i in myigenes if ( not isSignificant or (not(isnan(mypvals[i][isImputed])) and mypvals[i][isImputed] < 0.05) ) and ( not isImputed or (not isnan(myDATA[i][2])) ) ]

    print('isSignificant '+str(isSignificant)+', isImputed '+str(isImputed)+', '+geneSource+', '+area+', igenes = '+str(igenes)) 

    igenes_all.append(igenes[:])
    exprCoeffs_all.append([myDATA[i][1+isImputed] for i in igenes])
    igenes_in_exprData = [[iexprs[i] for i in range(0,len(foundGeneNames)) if foundGeneNames[i] == genes[igene]][0] for igene in igenes]
    exprData_all.append(array([exprData[i,:] for i in igenes_in_exprData]))
    isSignificants_all.append(isSignificant)
    isImputeds_all.append(isImputed)
    geneSources_all.append(geneSource)
    areas_all.append(area)


myVarStrICs = []
myVarValStrICs = []
myVarValStrsICs = []
myVarStrnonICs = []
myVarValStrnonICs = []
myVarValStrsnonICs = []
for igeneset in range(0,len(igenes_all)):
  print('isSignificant '+str(isSignificants_all[igeneset])+', isImputed '+str(isImputeds_all[igeneset])+', source='+geneSources_all[igeneset]+', area='+str(areas_all[igeneset]))
  print(str(list(zip([genes[i] for i in igenes_all[igeneset]],exprCoeffs_all[igeneset]))).replace('),','),\n'))
  varVals = []
  varVals_samps = []
  for ivar in range(0,len(variables)):
    varVals.append([])
    varVals_samps.append([])
  for igene in range(0,len(igenes_all[igeneset])):
    varVals[ivars_genes[igenes_all[igeneset][igene]]].append(exprCoeffs_all[igeneset][igene])
    varVals_samps[ivars_genes[igenes_all[igeneset][igene]]].append(exprData_all[igeneset][igene,:]/mean(exprData_all[igeneset][igene,isubjs_HC]))
  myVarStr = ''
  myVarValStr = ''
  myVarValStrs = []
  myVarStrIC = ''
  myVarValStrIC = ''
  myVarValStrsIC = [] 
  myVarStrnonIC = ''
  myVarValStrnonIC = ''
  myVarValStrsnonIC = []
  for i in range(0,exprData.shape[1]):
    myVarValStrs.append('')
    myVarValStrsIC.append('')
    myVarValStrsnonIC.append('')

  for ivar in range(0,len(variables)):
    if len(varVals[ivar]) > 0:
      myVarStr = myVarStr+variables[ivar]+','
      myVarValStr = myVarValStr+str(mean(varVals[ivar]))+','
      if variables[ivar] in ['gCa_HVAbar_Ca_HVA', 'gCa_LVAstbar_Ca_LVAst', 'gIhbar_Ih', 'gImbar_Im', 'gK_Pstbar_K_Pst', 'gNaT*_tbar_NaT*_t', 'gNap_Et2bar_Nap_Et2', 'g_pas']:
        if variables[ivar] == 'gNaT*_tbar_NaT*_t':
          myVarStrIC = myVarStrIC+'gNaTa_tbar_NaTa_t,gNaTs2_tbar_NaTs2_t,'
          myVarValStrIC = myVarValStrIC+'{:.3f}'.format(mean(varVals[ivar]))+','+'{:.3f}'.format(mean(varVals[ivar]))+','
          for isubj in range(0,exprData.shape[1]):
            myVarValStrsIC[isubj] = myVarValStrsIC[isubj]+'{:.3f}'.format(mean([varVals_samps[ivar][i][isubj] for i in range(0,len(varVals_samps[ivar]))]))+','+'{:.3f}'.format(mean([varVals_samps[ivar][i][isubj] for i in range(0,len(varVals_samps[ivar]))]))+','
        else:
          myVarStrIC = myVarStrIC+variables[ivar]+','
          myVarValStrIC = myVarValStrIC+'{:.3f}'.format(mean(varVals[ivar]))+','
          for isubj in range(0,exprData.shape[1]):
            myVarValStrsIC[isubj] = myVarValStrsIC[isubj]+'{:.3f}'.format(mean([varVals_samps[ivar][i][isubj] for i in range(0,len(varVals_samps[ivar]))]))+','
      else:
        if variables[ivar] == 'Calbin*':
          myVarStrnonIC = myVarStrnonIC+ 'Calbin,CalbinC,'
          myVarValStrnonIC = myVarValStrnonIC+'{:.3f}'.format(mean(varVals[ivar]))+','+'{:.3f}'.format(mean(varVals[ivar]))+','
          for isubj in range(0,exprData.shape[1]):
            myVarValStrsnonIC[isubj] = myVarValStrsnonIC[isubj]+'{:.3f}'.format(mean([varVals_samps[ivar][i][isubj] for i in range(0,len(varVals_samps[ivar]))]))+','+'{:.3f}'.format(mean([varVals_samps[ivar][i][isubj] for i in range(0,len(varVals_samps[ivar]))]))+','
        else:
          myVarStrnonIC = myVarStrnonIC+variables[ivar]+','
          myVarValStrnonIC = myVarValStrnonIC+'{:.3f}'.format(mean(varVals[ivar]))+','
          for isubj in range(0,exprData.shape[1]):
            myVarValStrsnonIC[isubj] = myVarValStrsnonIC[isubj]+'{:.3f}'.format(mean([varVals_samps[ivar][i][isubj] for i in range(0,len(varVals_samps[ivar]))]))+','
            
  print(myVarStr[0:-1]+' '+myVarValStr[0:-1])
  myVarStrICs.append(myVarStrIC)
  myVarValStrICs.append(myVarValStrIC)
  myVarValStrsICs.append(myVarValStrsIC[:])
  myVarStrnonICs.append(myVarStrnonIC)
  myVarValStrnonICs.append(myVarValStrnonIC)
  myVarValStrsnonICs.append(myVarValStrsnonIC[:])
  print(myVarStrnonIC[0:-1]+' '+myVarValStrnonIC[0:-1])
  print(myVarStrIC[0:-1]+' '+myVarValStrIC[0:-1])
  print(" ")

#for ivar in range(0,
for igeneset in range(0,len(myVarStrnonICs)):
  output_file = open('CMcomb_samps_nonIC_PFC_isSign'+str(isSignificants_all[igeneset])+'_isImputed'+str(isImputeds_all[igeneset])+'_'+geneSources_all[igeneset]+'.sh','w')
  output_file.write("""
BLOCKED='"""+myVarStrnonICs[igeneset][0:-1]+"""'\n
BLOCKEDCOEFFS=(""")
  for isubj in range(0,len(myVarValStrsnonICs[igeneset])):
    output_file.write(myVarValStrsnonICs[igeneset][isubj][0:-1]+'\n')
  output_file.write(')\n')

  output_file.close()

for igeneset in range(0,len(myVarStrICs)):
  output_file = open('CMcomb_samps_IC_PFC_isSign'+str(isSignificants_all[igeneset])+'_isImputed'+str(isImputeds_all[igeneset])+'_'+geneSources_all[igeneset]+'.sh','w')
  output_file.write("""
BLOCKED='"""+myVarStrICs[igeneset][0:-1]+"""'\n
BLOCKEDCOEFFS=(""")
  for isubj in range(0,len(myVarValStrsICs[igeneset])):
    output_file.write(myVarValStrsICs[igeneset][isubj][0:-1]+'\n')
  output_file.write(')\n')
  output_file.close()

#Checked that everything is correct:
#QWE=array([[float(x) for x in y[0:-1].split(',')] for y in myVarValStrsnonIC])
#[list(zip(genes,ACCDATA))[j] for j in igenes_all[5]]
#mean([QWE[i,:] for i in isubjs_SCZ],axis=0)/mean([QWE[i,:] for i in isubjs_HC],axis=0)

  
  
