
BLOCKED=GluR1,GluR1_memb,GluR2,GluR2_memb
BLOCKEDCOEFF=0.5,0.5,1.5,1.5
EXTRABLOCKED=Ca
EXTRABLOCKEDCOEFF=1.0
initfile=None

Econ=0.00015
wNMDA=1.0
gNap=0.03
Npulses=4
TRAINISIS=(-200.0 -180.0 -160.0 -140.0 -120.0 -100.0 -80.0 -60.0 -50.0 -40.0 -30.0 -25.0 -20.0 -15.0 -10.0 -5.0 0.0 5.0 10.0 15.0 20.0 25.0 30.0 40.0 50.0 60.0 80.0 100.0 120.0 140.0 160.0 180.0 200.0 -210.0 -220.0 -230.0 -2.5 2.5)
pulseamp=10.0
Nsyn=1
LOC=apic250-300

#Control simulations with and without neuromodulatory ligands:
for iISI in `seq 0 37`
do
  ISI=${TRAINISIS[iISI]}
  for LFLUX in 0.0 0.05
  do
    for ACHFLUX in 0.0 0.05
    do
      if [[ -f nrn_tstop25840000_${BLOCKED},${EXTRABLOCKED}x${BLOCKEDCOEFF},${EXTRABLOCKEDCOEFF}_n120_freq1.0_dur3.0_Lflux${LFLUX}_Gluflux0.0_AChflux${ACHFLUX}_Ntrains1_trainT100000.0_pair${ISI}_icell1_pulseamp10.0_Nsyn1_Econ0.00015_wNMDA1.0_gNap0.03_Npulses4_apic250-300.mat || -f nrn_tstop25840000_${BLOCKED},${EXTRABLOCKED}x${BLOCKEDCOEFF},${EXTRABLOCKEDCOEFF}_onset24040000.0_n120_freq1.0_dur3.0_Lflux${LFLUX}_Gluflux0.0_AChflux${ACHFLUX}_Ntrains1_trainT100000.0_pair${ISI}_icell1_pulseamp10.0_Nsyn1_Econ0.00015_wNMDA1.0_gNap0.03_Npulses4_apic250-300.mat ]]
      then
        echo "nrn_tstop25840000_${BLOCKED},${EXTRABLOCKED}x${BLOCKEDCOEFF},${EXTRABLOCKEDCOEFF}_n120_freq1.0_dur3.0_Lflux${LFLUX}_Gluflux0.0_AChflux${ACHFLUX}_Ntrains1_trainT100000.0_pair${ISI}_icell1_pulseamp10.0_Nsyn1_Econ0.00015_wNMDA1.0_gNap0.03_Npulses4_apic250-300.mat exists"
      else
        echo "nrn_tstop25840000_${BLOCKED},${EXTRABLOCKED}x${BLOCKEDCOEFF},${EXTRABLOCKEDCOEFF}_n120_freq1.0_dur3.0_Lflux${LFLUX}_Gluflux0.0_AChflux${ACHFLUX}_Ntrains1_trainT100000.0_pair${ISI}_icell1_pulseamp10.0_Nsyn1_Econ0.00015_wNMDA1.0_gNap0.03_Npulses4_apic250-300.mat does not exist"
        echo "python3 model_nrn_paired_withNap_contnm_var_npulses.py 25840000 1e-6 24040000 120 1.0 3 1.0 $LFLUX 0.0 $ACHFLUX 1 100000 $ISI 1 $Econ $wNMDA $gNap $Npulses $pulseamp $LOC $Nsyn None 23560000 600000 ${BLOCKED},${EXTRABLOCKED} ${BLOCKEDCOEFF},${EXTRABLOCKEDCOEFF}"
        python3 model_nrn_paired_withNap_contnm_var_npulses.py 25840000 1e-6 24040000 120 1.0 3 1.0 $LFLUX 0.0 $ACHFLUX 1 100000 $ISI 1 $Econ $wNMDA $gNap $Npulses $pulseamp $LOC $Nsyn None 23560000 600000 ${BLOCKED},${EXTRABLOCKED} ${BLOCKEDCOEFF},${EXTRABLOCKEDCOEFF}
      fi
    done
  done
done


#Simulations with synaptic CM variants, no ion-channel variants:
EXTRABLOCKEDS=('Calbin,CalbinC,DAGK,Gi,MGluR,NCX,Ng,PDE4,PKC,PP1' 'Calbin,CalbinC,DAGK,Gi,MGluR,NCX,Ng,PDE4,PKC,PP1' 'CK,DAGK,Gqabg,PKA' 'DAGK,PKA,PKC,PLA2,PP1' 'CK,Calbin,CalbinC,DAGK,Gi,Gqabg,MGluR,NCX,Ng,PDE4,PKA,PKC,PP1' 'Calbin,CalbinC,DAGK,Gi,MGluR,NCX,Ng,PDE4,PKA,PKC,PLA2,PP1' 'Calbin,CalbinC,DAGK,Gi,Ng,PDE4,PKC' 'DAGK,Gi,NCX,Ng,PDE4,PP1' 'CK,DAGK,Gqabg,PKA' 'DAGK,PKA,PKC,PLA2,PP1' 'CK,Calbin,CalbinC,DAGK,Gi,Gqabg,Ng,PDE4,PKA,PKC' 'DAGK,Gi,NCX,Ng,PDE4,PKA,PKC,PLA2,PP1' 'DAGK,MGluR,NCX,PDE4,PP1' 'MGluR' 'CK,DAGK,Gqabg,PKA' 'DAGK,PKA,PKC,PLA2,PP1' 'CK,DAGK,Gqabg,MGluR,NCX,PDE4,PKA,PP1' 'DAGK,MGluR,PKA,PKC,PLA2,PP1' 'DAGK,PDE4' 'CK,DAGK,Gqabg,PKA' 'DAGK,PKA,PKC,PLA2,PP1' 'CK,DAGK,Gqabg,PDE4,PKA' 'DAGK,PKA,PKC,PLA2,PP1')
EXTRABLOCKEDCOEFFS=(0.967,0.967,0.974,0.979,1.062,1.060,0.973,1.054,0.994,1.066 1.025,1.025,0.988,0.968,1.093,1.021,0.941,0.967,0.984,0.986 0.905,1.140,1.109,1.079 1.158,1.055,1.056,0.864,0.897 0.905,0.967,0.967,1.029,0.979,1.109,1.062,1.060,0.973,1.054,1.079,0.994,1.066 1.025,1.025,1.045,0.968,1.093,1.021,0.941,0.967,1.055,1.020,0.864,0.942 0.985,0.985,0.924,0.971,0.962,1.012,0.996 0.992,0.984,0.999,0.944,0.977,0.992 0.959,1.072,1.063,1.048 1.134,1.039,1.037,0.933,0.940 0.959,0.985,0.985,0.998,0.971,1.063,0.962,1.012,1.048,0.996 1.039,0.984,0.999,0.944,0.977,1.039,1.037,0.933,0.966 0.918,1.062,1.090,1.054,1.066 1.093 0.905,1.140,1.109,1.079 1.158,1.055,1.056,0.864,0.897 0.905,1.029,1.109,1.062,1.090,1.054,1.079,1.066 1.158,1.093,1.055,1.056,0.864,0.897 0.924,1.012 0.959,1.072,1.063,1.048 1.134,1.039,1.037,0.933,0.940 0.959,0.998,1.063,1.012,1.048 1.134,1.039,1.037,0.933,0.940)
imut=0 #no ion-channel variants
for iISI in `seq 0 37`
do
  ISI=${TRAINISIS[iISI]}
  for iextra in 10 11
  do
    EXTRABLOCKED=${EXTRABLOCKEDS[iextra]}
    EXTRABLOCKEDCOEFF=${EXTRABLOCKEDCOEFFS[iextra]}
    extfilename=nrn_paired_imutCMcomb-1_comb${iextra}_${Econ}_${wNMDA}_${ISI}.mat
    if [ -f $extfilename ]
    then
      echo "$extfilename exists"
    else
      echo "$extfilename does not exist"
      echo "python3 model_nrn_paired_withNap_contnm_var_npulses_muts_extfilename.py 25840000 1e-6 24040000 120 1.0 3 1.0 0.05 0.0 0.05 1 100000 $ISI 1 $Econ $wNMDA $gNap $imut $Npulses $pulseamp $LOC $Nsyn None 23560000 600000 ${BLOCKED},${EXTRABLOCKED} ${BLOCKEDCOEFF},${EXTRABLOCKEDCOEFF} $ALTERED $ALTEREDCOEFF $extfilename"
      python3 model_nrn_paired_withNap_contnm_var_npulses_muts_extfilename.py 25840000 1e-6 24040000 120 1.0 3 1.0 0.05 0.0 0.05 1 100000 $ISI 1 $Econ $wNMDA $gNap $imut $Npulses $pulseamp $LOC $Nsyn None 23560000 600000 ${BLOCKED},${EXTRABLOCKED} ${BLOCKEDCOEFF},${EXTRABLOCKEDCOEFF} $ALTERED $ALTEREDCOEFF $extfilename
    fi
  done
done

#Simulations with ion-channel variants, no synaptic variants
EXTRABLOCKED=Ca
EXTRABLOCKEDCOEFF=1.0
for iISI in `seq 0 37`
do
  ISI=${TRAINISIS[iISI]}
  for imut in 17 20 23 32 35 38 41 44 47
  do
    if [ -f nrn_tstop25840000_${BLOCKED},${EXTRABLOCKED}x${BLOCKEDCOEFF},${EXTRABLOCKEDCOEFF}_n120_freq1.0_dur3.0_Lflux0.05_Gluflux0.0_AChflux0.05_Ntrains1_trainT100000.0_pair${ISI}_icell1_imutc${imut}_pulseamp10.0_Nsyn1_Econ0.00015_wNMDA1.0_gNap0.03_Npulses4_apic250-300.mat ]
    then
      echo "nrn_tstop25840000_${BLOCKED},${EXTRABLOCKED}x${BLOCKEDCOEFF},${EXTRABLOCKEDCOEFF}_n120_freq1.0_dur3.0_Lflux0.05_Gluflux0.0_AChflux0.05_Ntrains1_trainT100000.0_pair${ISI}_icell1_imutc${imut}_pulseamp10.0_Nsyn1_Econ0.00015_wNMDA1.0_gNap0.03_Npulses4_apic250-300.mat exists"
    else
      echo "nrn_tstop25840000_${BLOCKED},${EXTRABLOCKED}x${BLOCKEDCOEFF},${EXTRABLOCKEDCOEFF}_n120_freq1.0_dur3.0_Lflux0.05_Gluflux0.0_AChflux0.05_Ntrains1_trainT100000.0_pair${ISI}_icell1_imutc${imut}_pulseamp10.0_Nsyn1_Econ0.00015_wNMDA1.0_gNap0.03_Npulses4_apic250-300.mat does not exist"
      echo "python3 model_nrn_paired_withNap_contnm_var_npulses_muts_extfilename.py 25840000 1e-6 24040000 120 1.0 3 1.0 0.05 0.0 0.05 1 100000 $ISI 1 $Econ $wNMDA $gNap $imut $Npulses $pulseamp $LOC $Nsyn None 23560000 600000 ${BLOCKED},${EXTRABLOCKED} ${BLOCKEDCOEFF},${EXTRABLOCKEDCOEFF}"
      python3 model_nrn_paired_withNap_contnm_var_npulses_muts_extfilename.py 25840000 1e-6 24040000 120 1.0 3 1.0 0.05 0.0 0.05 1 100000 $ISI 1 $Econ $wNMDA $gNap $imut $Npulses $pulseamp $LOC $Nsyn None 23560000 600000 ${BLOCKED},${EXTRABLOCKED} ${BLOCKEDCOEFF},${EXTRABLOCKEDCOEFF}
    fi
  done
done

#Simulations with ion-channel CM variants, no synaptic variants
for iISI in `seq 0 37`
do
  ISI=${TRAINISIS[iISI]}
  for imut in 10 11
  do
    extfilename=nrn_paired_imutCMcomb${imut}_comb-1_${Econ}_${wNMDA}_${ISI}.mat
    if [ -f $extfilename ]
    then
      echo "$extfilename exists"
    else
      echo "$extfilename does not exist"
      echo "python3 model_nrn_paired_withNap_contnm_var_npulses_mutCMcombs_extfilename.py 25840000 1e-6 24040000 120 1.0 3 1.0 0.05 0.0 0.05 1 100000 $ISI 1 $Econ $wNMDA $gNap $imut $Npulses $pulseamp $LOC $Nsyn None 23560000 600000 ${BLOCKED} ${BLOCKEDCOEFF} $ALTERED $ALTEREDCOEFF $extfilename"
      python3 model_nrn_paired_withNap_contnm_var_npulses_mutCMcombs_extfilename.py 25840000 1e-6 24040000 120 1.0 3 1.0 0.05 0.0 0.05 1 100000 $ISI 1 $Econ $wNMDA $gNap $imut $Npulses $pulseamp $LOC $Nsyn None 23560000 600000 ${BLOCKED} ${BLOCKEDCOEFF} $ALTERED $ALTEREDCOEFF $extfilename
    fi
  done
done

#Simulations with ion-channel CM variants, with synaptic CM variants
for iISI in `seq 0 37`
do
  ISI=${TRAINISIS[iISI]}
  for imut in 10 11
  do
    extfilename=nrn_paired_imutCMcomb${imut}_comb${imut}_${Econ}_${wNMDA}_${ISI}.mat
    if [ -f $extfilename ]
    then
      echo "$extfilename exists"
    else
      echo "$extfilename does not exist"
      echo "python3 model_nrn_paired_withNap_contnm_var_npulses_mutCMcombs_extfilename.py 25840000 1e-6 24040000 120 1.0 3 1.0 0.05 0.0 0.05 1 100000 $ISI 1 $Econ $wNMDA $gNap $imut $Npulses $pulseamp $LOC $Nsyn None 23560000 600000 ${BLOCKED},${EXTRABLOCKED} ${BLOCKEDCOEFF},${EXTRABLOCKEDCOEFF} $ALTERED $ALTEREDCOEFF $extfilename"
      python3 model_nrn_paired_withNap_contnm_var_npulses_mutCMcombs_extfilename.py 25840000 1e-6 24040000 120 1.0 3 1.0 0.05 0.0 0.05 1 100000 $ISI 1 $Econ $wNMDA $gNap $imut $Npulses $pulseamp $LOC $Nsyn None 23560000 600000 ${BLOCKED},${EXTRABLOCKED} ${BLOCKEDCOEFF},${EXTRABLOCKEDCOEFF} $ALTERED $ALTEREDCOEFF $extfilename
    fi
  done
done
