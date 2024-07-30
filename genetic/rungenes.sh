module load plink/1.90b6.2
for iset in 0 1 2 3 4 5
do
setNames=(genes_new genes_new_synaptic genes_all_synaptic_PKA genes_all_synaptic_PKC genes_new_ionchannels genes_lips_synaptic_ABFGJLOQ)
pwd
for THR in 1e-06 5e-06 1e-05
do
  myset=${setNames[iset]}
  cd norment_i/bed
  plink --bfile all --extract ~/analysis/snps_of_interest_${myset}_${THR}.txt --out ~/analysis/genetic/${myset}_${THR}_i --keep ~/analysis/subj_IDs_of_interest.txt --recode
  cd norment_ii/bed
  plink --bfile all --extract ~/analysis/snps_of_interest_${myset}_${THR}.txt --out ~/analysis/genetic/${myset}_${THR}_ii --keep ~/analysis/subj_IDs_of_interest.txt --recode
  cd norment_iii/bed
  plink --bfile all --extract ~/analysis/snps_of_interest_${myset}_${THR}.txt --out ~/analysis/genetic/${myset}_${THR}_iii --keep ~/analysis/subj_IDs_of_interest.txt --recode

  cd ~/analysis/genetic
  for q in ${myset}_${THR}_i ${myset}_${THR}_ii ${myset}_${THR}_iii
  do
      echo "plink --file ${q} --make-bed --out ${q}_new"
      plink --file ${q} --make-bed --out ${q}_new
  done
  
  
  for q in ${myset}_${THR}
  do
    ls ${q}_i*_new*bed |sed 's/.bed//' > mergelist${myset}.txt
    echo "plink --merge-list mergelist${myset}.txt --make-bed --out ${q}_merged"
    plink --merge-list mergelist${myset}.txt --make-bed --out ${q}_merged
  done

  ~/PRSice_linux --base ~/PGC_SCZ_0518_EUR_noTOP.sumstats --beta --target ${myset}_${THR}_merged --pheno DUMMY_PHENO.txt --stat BETA --pvalue PVAL --print-snp --no-full --snp VARIANT_ID --bar-levels "1" --fastscore --out ${myset}_${THR}_PRS --no-clump --maf 0.01 --thread 8
  ~/PRSice_linux --base ~/PGC_SCZ_0518_EUR_noTOP.sumstats --beta --target ${myset}_${THR}_merged --pheno DUMMY_PHENO.txt --stat BETA --pvalue PVAL --print-snp --no-full --snp VARIANT_ID --bar-levels "1" --fastscore --out ${myset}_${THR}_PRS --no-clump --maf 0.01 --thread 8 --extract ${myset}_${THR}_PRS.valid
  echo "rm ${myset}_${THR}_PRS.valid"
  rm ${myset}_${THR}_PRS.valid

  cat ${myset}_${THR}_PRS.best |cut -f1 -d' ' > subjs_${myset}_${THR}.txt 
done

done
