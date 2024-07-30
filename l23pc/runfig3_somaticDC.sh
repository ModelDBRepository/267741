#!/bin/bash
# Job name:
#SBATCH --job-name=l23pc
#
# Project:
#SBATCH --account=nn9529k
#
# Wall clock limit:
#SBATCH --time=36:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=1600M
#SBATCH --nodes=1 --ntasks-per-node=1

module load matplotlib/3.4.2-foss-2021a

ICELL=0
ATOL=0.00005

AMPS=(0.0 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25)

for iamp in `seq 0 10`
do
  IMUT=0
  echo "python3 runmodelmut_somaticDC.py $ICELL $IMUT ${AMPS[iamp]} $ATOL"
  python3 runmodelmut_somaticDC.py $ICELL $IMUT ${AMPS[iamp]} $ATOL

  for IMUT in 2 5 8 11 14 17 20 23 26 29 32 35 38 41 44 47
  do
    echo "python3 runmodelmutc_somaticDC.py $ICELL $IMUT ${AMPS[iamp]} $ATOL"
    python3 runmodelmutc_somaticDC.py $ICELL $IMUT ${AMPS[iamp]} $ATOL
  done

  for IMUT in 10 11
  do
    echo "python3 runmodelmutCMcomb_somaticDC.py $ICELL $IMUT ${AMPS[iamp]} $ATOL"
    python3 runmodelmutCMcomb_somaticDC.py $ICELL $IMUT ${AMPS[iamp]} $ATOL
  done
    
done

