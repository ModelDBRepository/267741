#!/bin/bash
# Job name:
#SBATCH --job-name=l23pc
#
# Project:
#SBATCH --account=nn9529k
#
# Wall clock limit:
#SBATCH --time=48:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=3200M
#SBATCH --nodes=1 --ntasks-per-node=1

module load matplotlib/3.4.2-foss-2021a

ICELL=0
RATES=(0.7)
FREQ=1.0
NSTIM=100
ISIS=(-120.0 -100.0 -80.0 -60.0 -50.0 -40.0 -30.0 -25.0 -20.0 -15.0 -10.0 -5.0 -2.5 0.0 2.5 5.0 10.0 15.0 20.0 25.0 30.0 40.0 50.0 60.0 80.0 100.0 120.0 140.0)

#PULSEAMPS=(0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7)
NINPUTS=1
ECON=0.00015
WNMDA=1.0
NSAMP=200
Npulses=4

DENDTREE=apic
SPINELOCATIONS=250-300
NSYN=1
RATE=0.7

neckLen=0.5
neckDiam=0.1
NSYN=1
PULSEAMP=10.0
gNap=0.03

for iISI in `seq 0 27`
do
  #These simulations took me around 12 hours per ISI value. Parallelization recommended.
  ISI=${ISIS[iISI]}

  #1) Simulations without alterations of ion-channel conductances
  IMUT=0
  if [ -f currClips${ICELL}_imut${IMUT}_neckLen${neckLen}_neckDiam${neckDiam}_stimfreq${FREQ}_pulseamp${PULSEAMP}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_Npulses${Npulses}_ISI${ISI}.mat ]
  then
    echo "currClips${ICELL}_imut${IMUT}_neckLen${neckLen}_neckDiam${neckDiam}_stimfreq${FREQ}_pulseamp${PULSEAMP}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_Npulses${Npulses}_ISI${ISI}.mat already exists"
  else
    echo "currClips${ICELL}_imut${IMUT}_neckLen${neckLen}_neckDiam${neckDiam}_stimfreq${FREQ}_pulseamp${PULSEAMP}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_Npulses${Npulses}_ISI${ISI}.mat does not exist"
    for RDSEED in `seq 0 $((NSAMP-1))`
    do
      if [ -f noisy_icell${ICELL}_imut${IMUT}_n${NSTIM}_${FREQ}_neckLen${neckLen}_neckDiam${neckDiam}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_rateE${RATE}_Npulses${Npulses}_ISI${ISI}_dtpulses10.0_pulseamp${PULSEAMP}_seed${RDSEED}.mat ]
      then
        echo "noisy_icell${ICELL}_imut${IMUT}_n${NSTIM}_${FREQ}_neckLen${neckLen}_neckDiam${neckDiam}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_rateE${RATE}_Npulses${Npulses}_ISI${ISI}_dtpulses10.0_pulseamp${PULSEAMP}_seed${RDSEED}.mat  exists"
      else
        echo "noisy_icell${ICELL}_imut${IMUT}_n${NSTIM}_${FREQ}_neckLen${neckLen}_neckDiam${neckDiam}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_rateE${RATE}_Npulses${Npulses}_ISI${ISI}_dtpulses10.0_pulseamp${PULSEAMP}_seed${RDSEED}.mat does not exist"
        echo "python3 runmodelmut_withNap.py $ICELL $IMUT ${FREQ} ${NSTIM} $NSYN $NINPUTS ${DENDTREE} ${SPINELOCATIONS} $RATE $Npulses $ISI 10.0 ${PULSEAMP} ${neckLen} ${neckDiam} $ECON $WNMDA $gNap $RDSEED 0"
        python3 runmodelmut_withNap.py $ICELL $IMUT ${FREQ} ${NSTIM} $NSYN $NINPUTS ${DENDTREE} ${SPINELOCATIONS} $RATE $Npulses $ISI 10.0 ${PULSEAMP} ${neckLen} ${neckDiam} $ECON $WNMDA $gNap $RDSEED 0
      fi
    done

    sleep 6

    echo "python3 collectlocalcurrsmutonepulse_withNap.py $ICELL $IMUT $NSYN $ECON $WNMDA $gNap $PULSEAMP $DENDTREE $SPINELOCATIONS ${FREQ} $ISI $Npulses $neckLen $neckDiam"
    python3 collectlocalcurrsmutonepulse_withNap.py $ICELL $IMUT $NSYN $ECON $WNMDA $gNap $PULSEAMP $DENDTREE $SPINELOCATIONS ${FREQ} $ISI $Npulses $neckLen $neckDiam

    if [ -f currClips${ICELL}_imut${IMUT}_neckLen${neckLen}_neckDiam${neckDiam}_stimfreq${FREQ}_pulseamp${PULSEAMP}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_Npulses${Npulses}_ISI${ISI}.mat ]
    then
      echo "rm noisy_icell${ICELL}_imut${IMUT}_n${NSTIM}_${FREQ}_neckLen${neckLen}_neckDiam${neckDiam}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_rateE${RATE}_Npulses${Npulses}_ISI${ISI}_dtpulses10.0_pulseamp${PULSEAMP}_seed*.mat"
      rm noisy_icell${ICELL}_imut${IMUT}_n${NSTIM}_${FREQ}_neckLen${neckLen}_neckDiam${neckDiam}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_rateE${RATE}_Npulses${Npulses}_ISI${ISI}_dtpulses10.0_pulseamp${PULSEAMP}_seed*.mat
    else
      echo "File currClips${ICELL}_imut${IMUT}_neckLen${neckLen}_neckDiam${neckDiam}_stimfreq${FREQ}_pulseamp${PULSEAMP}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_Npulses${Npulses}_ISI${ISI}.mat not created"
      echo "Something wrong, currClips not saved"
    fi
  fi

  #2) Simulations with alterations of single ion-channel conductances (+-20%)
  for IMUT in 2 5 8 11 14 17 20 23 26 29 32 35 38 41 44 47
  do	      
    if [ -f currClips${ICELL}_imutc${IMUT}_neckLen${neckLen}_neckDiam${neckDiam}_stimfreq${FREQ}_pulseamp${PULSEAMP}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_Npulses${Npulses}_ISI${ISI}.mat ]
    then
      echo "currClips${ICELL}_imutc${IMUT}_neckLen${neckLen}_neckDiam${neckDiam}_stimfreq${FREQ}_pulseamp${PULSEAMP}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_Npulses${Npulses}_ISI${ISI}.mat already exists"
    else
      echo "currClips${ICELL}_imutc${IMUT}_neckLen${neckLen}_neckDiam${neckDiam}_stimfreq${FREQ}_pulseamp${PULSEAMP}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_Npulses${Npulses}_ISI${ISI}.mat does not exist"
      for RDSEED in `seq 0 $((NSAMP-1))`
      do
        if [ -f noisy_icell${ICELL}_imutc${IMUT}_n${NSTIM}_${FREQ}_neckLen${neckLen}_neckDiam${neckDiam}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_rateE${RATE}_Npulses${Npulses}_ISI${ISI}_dtpulses10.0_pulseamp${PULSEAMP}_seed${RDSEED}.mat ]
        then
          echo "noisy_icell${ICELL}_imutc${IMUT}_n${NSTIM}_${FREQ}_neckLen${neckLen}_neckDiam${neckDiam}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_rateE${RATE}_Npulses${Npulses}_ISI${ISI}_dtpulses10.0_pulseamp${PULSEAMP}_seed${RDSEED}.mat  exists"
        else
          echo "noisy_icell${ICELL}_imutc${IMUT}_n${NSTIM}_${FREQ}_neckLen${neckLen}_neckDiam${neckDiam}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_rateE${RATE}_Npulses${Npulses}_ISI${ISI}_dtpulses10.0_pulseamp${PULSEAMP}_seed${RDSEED}.mat does not exist"
          echo "python3 runmodelmutc_withNap.py $ICELL $IMUT ${FREQ} ${NSTIM} $NSYN $NINPUTS ${DENDTREE} ${SPINELOCATIONS} $RATE $Npulses $ISI 10.0 ${PULSEAMP} ${neckLen} ${neckDiam} $ECON $WNMDA $gNap $RDSEED 0"
          python3 runmodelmutc_withNap.py $ICELL $IMUT ${FREQ} ${NSTIM} $NSYN $NINPUTS ${DENDTREE} ${SPINELOCATIONS} $RATE $Npulses $ISI 10.0 ${PULSEAMP} ${neckLen} ${neckDiam} $ECON $WNMDA $gNap $RDSEED 0
        fi
      done

      sleep 6

      echo "python3 collectlocalcurrsmutconepulse_withNap.py $ICELL $IMUT $NSYN $ECON $WNMDA $gNap $PULSEAMP $DENDTREE $SPINELOCATIONS ${FREQ} $ISI $Npulses $neckLen $neckDiam"
      python3 collectlocalcurrsmutconepulse_withNap.py $ICELL $IMUT $NSYN $ECON $WNMDA $gNap $PULSEAMP $DENDTREE $SPINELOCATIONS ${FREQ} $ISI $Npulses $neckLen $neckDiam

      if [ -f currClips${ICELL}_imutc${IMUT}_neckLen${neckLen}_neckDiam${neckDiam}_stimfreq${FREQ}_pulseamp${PULSEAMP}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_Npulses${Npulses}_ISI${ISI}.mat ]
      then
        echo "rm noisy_icell${ICELL}_imutc${IMUT}_n${NSTIM}_${FREQ}_neckLen${neckLen}_neckDiam${neckDiam}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_rateE${RATE}_Npulses${Npulses}_ISI${ISI}_dtpulses10.0_pulseamp${PULSEAMP}_seed*.mat"
        rm noisy_icell${ICELL}_imutc${IMUT}_n${NSTIM}_${FREQ}_neckLen${neckLen}_neckDiam${neckDiam}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_rateE${RATE}_Npulses${Npulses}_ISI${ISI}_dtpulses10.0_pulseamp${PULSEAMP}_seed*.mat
      else
        echo "File currClips${ICELL}_imutc${IMUT}_neckLen${neckLen}_neckDiam${neckDiam}_stimfreq${FREQ}_pulseamp${PULSEAMP}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_Npulses${Npulses}_ISI${ISI}.mat not created"
        echo "Something wrong, currClips not saved"
      fi
    fi
  done
  
  #3) Simulations with combinations of alterations of ion-channel conductances
  for IMUT in 10 11
  do
    if [ -f currClips${ICELL}_imutCMcomb${IMUT}_neckLen${neckLen}_neckDiam${neckDiam}_stimfreq${FREQ}_pulseamp${PULSEAMP}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_Npulses${Npulses}_ISI${ISI}.mat ]
    then
      echo "currClips${ICELL}_imutCMcomb${IMUT}_neckLen${neckLen}_neckDiam${neckDiam}_stimfreq${FREQ}_pulseamp${PULSEAMP}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_Npulses${Npulses}_ISI${ISI}.mat already exists"
    else
      echo "currClips${ICELL}_imutCMcomb${IMUT}_neckLen${neckLen}_neckDiam${neckDiam}_stimfreq${FREQ}_pulseamp${PULSEAMP}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_Npulses${Npulses}_ISI${ISI}.mat does not exist"
      for RDSEED in `seq 0 $((NSAMP-1))`
      do
        if [ -f noisy_icell${ICELL}_imutCMcomb${IMUT}_n${NSTIM}_${FREQ}_neckLen${neckLen}_neckDiam${neckDiam}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_rateE${RATE}_Npulses${Npulses}_ISI${ISI}_dtpulses10.0_pulseamp${PULSEAMP}_seed${RDSEED}.mat ]
        then
          echo "noisy_icell${ICELL}_imutCMcomb${IMUT}_n${NSTIM}_${FREQ}_neckLen${neckLen}_neckDiam${neckDiam}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_rateE${RATE}_Npulses${Npulses}_ISI${ISI}_dtpulses10.0_pulseamp${PULSEAMP}_seed${RDSEED}.mat  exists"
        else
          echo "noisy_icell${ICELL}_imutCMcomb${IMUT}_n${NSTIM}_${FREQ}_neckLen${neckLen}_neckDiam${neckDiam}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_rateE${RATE}_Npulses${Npulses}_ISI${ISI}_dtpulses10.0_pulseamp${PULSEAMP}_seed${RDSEED}.mat does not exist"
	  echo "python3 runmodelmutCMcomb_withNap.py $ICELL $IMUT ${FREQ} ${NSTIM} $NSYN $NINPUTS ${DENDTREE} ${SPINELOCATIONS} $RATE $Npulses $ISI 10.0 ${PULSEAMP} ${neckLen} ${neckDiam} $ECON $WNMDA $gNap $RDSEED 0"
	  python3 runmodelmutCMcomb_withNap.py $ICELL $IMUT ${FREQ} ${NSTIM} $NSYN $NINPUTS ${DENDTREE} ${SPINELOCATIONS} $RATE $Npulses $ISI 10.0 ${PULSEAMP} ${neckLen} ${neckDiam} $ECON $WNMDA $gNap $RDSEED 0
        fi
      done

      sleep 6

      echo "python3 collectlocalcurrsmutCMcombonepulse_withNap.py $ICELL $IMUT $NSYN $ECON $WNMDA $gNap $PULSEAMP $DENDTREE $SPINELOCATIONS ${FREQ} $ISI $Npulses $neckLen $neckDiam"
      python3 collectlocalcurrsmutCMcombonepulse_withNap.py $ICELL $IMUT $NSYN $ECON $WNMDA $gNap $PULSEAMP $DENDTREE $SPINELOCATIONS ${FREQ} $ISI $Npulses $neckLen $neckDiam

      if [ -f currClips${ICELL}_imutCMcomb${IMUT}_neckLen${neckLen}_neckDiam${neckDiam}_stimfreq${FREQ}_pulseamp${PULSEAMP}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_Npulses${Npulses}_ISI${ISI}.mat ]
      then
        echo "rm noisy_icell${ICELL}_imutCMcomb${IMUT}_n${NSTIM}_${FREQ}_neckLen${neckLen}_neckDiam${neckDiam}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_rateE${RATE}_Npulses${Npulses}_ISI${ISI}_dtpulses10.0_pulseamp${PULSEAMP}_seed*.mat"
	rm noisy_icell${ICELL}_imutCMcomb${IMUT}_n${NSTIM}_${FREQ}_neckLen${neckLen}_neckDiam${neckDiam}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_rateE${RATE}_Npulses${Npulses}_ISI${ISI}_dtpulses10.0_pulseamp${PULSEAMP}_seed*.mat
      else
        echo "File currClips${ICELL}_imutCMcomb${IMUT}_neckLen${neckLen}_neckDiam${neckDiam}_stimfreq${FREQ}_pulseamp${PULSEAMP}_Nsyn${NSYN}_Ninputs${NINPUTS}${DENDTREE}${SPINELOCATIONS}_Econ${ECON}_wNMDA${WNMDA}_gNap${gNap}_Npulses${Npulses}_ISI${ISI}.mat not created"
	echo "Something wrong, currClips not saved"
      fi
    fi
  done
done

