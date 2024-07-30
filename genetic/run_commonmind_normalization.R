##Instructions
#Postmortem gene expression data from prefrontal cortex are available through the CommonMind Knowledge Portal (https://www.synapse.org/#!Synapse:syn2759792/wiki/197295).
#The raw data from the ACC is found here: https://www.synapse.org/#!Synapse:syn29442240
#The raw data from the PFC is found here: https://www.synapse.org/#!Synapse:syn18097440
#You will need to sign a distribution agreement to access these data. Once you have accessed and downloaded the data, collect them into a countData matrix, after which the following script can be run (this should be done separately for ACC and PFC data).

## pre-filter
library(edgeR)
cpm <- edgeR::cpm(countData)
l <- nrow(countData)*0.5
keep <- rowSums(cpm>1)>=l # require more than 1 cpm in at least 50% of samples
data.filt <- data[keep,]

## normalization
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=data.filt, colData=metaData, design= ~ Diagnosis) 
dds <- DESeq(dds)
normData <- counts(dds, normalized=TRUE)

##Post-processing
#The resulting data should be saved and imputed to get estimates of neuronal RNA expression in each subject. The imputation was done using CIBERSORTx with the following parameters:
#-Reference file: Zhang et al. 2016 "Purification and Characterization of Progenitor and Mature Human Astrocytes Reveals Transcriptional and Functional Differences with Mouse", Neuron 6;89(1):37-53, Supplementary material 3, Human only
#-Window size 20
#-High resolution
#-Batch correction enabled
#The imputed data should be saved to CIBERSORTxHiRes_ACC_Neurons_Window20.txt and CIBERSORTxHiRes_PFC_Neurons_Window20.txt, these will be further processed by extract_commonmind_data_ACC.py and extract_commonmind_data_PFC.py


