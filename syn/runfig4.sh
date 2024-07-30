
GLUR=GluR1,GluR1_memb,GluR2,GluR2_memb
ALTERED=1
ALTEREDCOEFF=1.0

BLOCKEDS=('Gs' 'MGluR' 'PP1' 'PDE4' 'PP2A' 'DAGK' 'Gi' 'Calbin' 'PKC' 'PLC' 'NCX' 'PP2B' 'PKA' 'Ng' 'Gqabg' 'CK' 'PLA2' 'I1')
BLOCKEDCOEFFS=(0.8 1.2)
BLOCKEDS2=('Calbin,CalbinC')
BLOCKEDCOEFFS2=(0.8,0.8 1.2,1.2)
BLOCKEDCOMBS=('Calbin,CalbinC,DAGK,Gi,MGluR,NCX,Ng,PDE4,PKC,PP1' 'Calbin,CalbinC,DAGK,Gi,MGluR,NCX,Ng,PDE4,PKC,PP1' 'CK,DAGK,Gqabg,PKA' 'DAGK,PKA,PKC,PLA2,PP1' 'CK,Calbin,CalbinC,DAGK,Gi,Gqabg,MGluR,NCX,Ng,PDE4,PKA,PKC,PP1' 'Calbin,CalbinC,DAGK,Gi,MGluR,NCX,Ng,PDE4,PKA,PKC,PLA2,PP1' 'Calbin,CalbinC,DAGK,Gi,Ng,PDE4,PKC' 'DAGK,Gi,NCX,Ng,PDE4,PP1' 'CK,DAGK,Gqabg,PKA' 'DAGK,PKA,PKC,PLA2,PP1' 'CK,Calbin,CalbinC,DAGK,Gi,Gqabg,Ng,PDE4,PKA,PKC' 'DAGK,Gi,NCX,Ng,PDE4,PKA,PKC,PLA2,PP1' 'DAGK,MGluR,NCX,PDE4,PP1' 'MGluR' 'CK,DAGK,Gqabg,PKA' 'DAGK,PKA,PKC,PLA2,PP1' 'CK,DAGK,Gqabg,MGluR,NCX,PDE4,PKA,PP1' 'DAGK,MGluR,PKA,PKC,PLA2,PP1' 'DAGK,PDE4' 'CK,DAGK,Gqabg,PKA' 'DAGK,PKA,PKC,PLA2,PP1' 'CK,DAGK,Gqabg,PDE4,PKA' 'DAGK,PKA,PKC,PLA2,PP1')
BLOCKEDCOMBCOEFFS=(0.967,0.967,0.974,0.979,1.062,1.060,0.973,1.054,0.994,1.066 1.025,1.025,0.988,0.968,1.093,1.021,0.941,0.967,0.984,0.986 0.905,1.140,1.109,1.079 1.158,1.055,1.056,0.864,0.897 0.905,0.967,0.967,1.029,0.979,1.109,1.062,1.060,0.973,1.054,1.079,0.994,1.066 1.025,1.025,1.045,0.968,1.093,1.021,0.941,0.967,1.055,1.020,0.864,0.942 0.985,0.985,0.924,0.971,0.962,1.012,0.996 0.992,0.984,0.999,0.944,0.977,0.992 0.959,1.072,1.063,1.048 1.134,1.039,1.037,0.933,0.940 0.959,0.985,0.985,0.998,0.971,1.063,0.962,1.012,1.048,0.996 1.039,0.984,0.999,0.944,0.977,1.039,1.037,0.933,0.966 0.918,1.062,1.090,1.054,1.066 1.093 0.905,1.140,1.109,1.079 1.158,1.055,1.056,0.864,0.897 0.905,1.029,1.109,1.062,1.090,1.054,1.079,1.066 1.158,1.093,1.055,1.056,0.864,0.897 0.924,1.012 0.959,1.072,1.063,1.048 1.134,1.039,1.037,0.933,0.940 0.959,0.998,1.063,1.012,1.048 1.134,1.039,1.037,0.933,0.940)

initfile=None

TSHORT=32000000
ONSET=24040000

CAFLUX=150.0
LFLUX=5.0
GLUFLUX=10.0
ACHFLUX=10.0
initfile=None

Ntrains=1200
for GLURCOEFF in 2.0,2.0,0.0,0.0 0.0,0.0,2.0,2.0
do
  #1) Siimulations without NM
  if [ -f nrn_tstop32000000_tol1e-06_${GLUR}x${GLURCOEFF}_k1x1.0_onset24040000.0_n1_freq100.0_dur3.0_flux${CAFLUX}_Lflux0.0_Gluflux0.0_AChflux0.0_Ntrains1200_trainT500.0.mat ]
  then
    echo "nrn_tstop32000000_tol1e-06_${GLUR}x${GLURCOEFF}_k1x1.0_onset24040000.0_n1_freq100.0_dur3.0_flux${CAFLUX}_Lflux0.0_Gluflux0.0_AChflux0.0_Ntrains1200_trainT500.0.mat exists"
  else
    echo "nrn_tstop32000000_tol1e-06_${GLUR}x${GLURCOEFF}_k1x1.0_onset24040000.0_n1_freq100.0_dur3.0_flux${CAFLUX}_Lflux0.0_Gluflux0.0_AChflux0.0_Ntrains1200_trainT500.0.mat does not exist"
    echo "python3 model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 1 100 3.0 $CAFLUX 0.0 0.0 0.0 $Ntrains 500 $initfile $GLUR $GLURCOEFF $ALTERED $ALTEREDCOEFF"
    python3 model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 1 100 3.0 $CAFLUX 0.0 0.0 0.0 $Ntrains 500 $initfile $GLUR $GLURCOEFF $ALTERED $ALTEREDCOEFF
  fi

  #2) Simulations with NM
  if [ -f nrn_tstop32000000_tol1e-06_${GLUR}x${GLURCOEFF}_k1x1.0_onset24040000.0_n1_freq100.0_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1200_trainT500.0.mat ]
  then
    echo "nrn_tstop32000000_tol1e-06_${GLUR}x${GLURCOEFF}_k1x1.0_onset24040000.0_n1_freq100.0_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1200_trainT500.0.mat exists"
  else
    echo "nrn_tstop32000000_tol1e-06_${GLUR}x${GLURCOEFF}_k1x1.0_onset24040000.0_n1_freq100.0_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1200_trainT500.0.mat does not exist"
    echo "python3 model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 1 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX $Ntrains 500 $initfile $GLUR $GLURCOEFF $ALTERED $ALTEREDCOEFF"
    python3 model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 1 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX $Ntrains 500 $initfile $GLUR $GLURCOEFF $ALTERED $ALTEREDCOEFF
  fi
  
  #3) Simulations with single protein concentration changed
  for iblockcoeff in 0 1
  do
    for iblock in `seq 0 16`
    do
      BLOCKED=${GLUR},${BLOCKEDS[iblock]}
      BLOCKEDCOEFF=${GLURCOEFF},${BLOCKEDCOEFFS[iblockcoeff]}      
      if [ -f nrn_tstop32000000_tol1e-06_${BLOCKED}x${BLOCKEDCOEFF}_k1x1.0_onset24040000.0_n1_freq100.0_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1200_trainT500.0.mat ]
      then
        echo "nrn_tstop32000000_tol1e-06_${BLOCKED}x${BLOCKEDCOEFF}_k1x1.0_onset24040000.0_n1_freq100.0_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1200_trainT500.0.mat exists"
      else
        echo "nrn_tstop32000000_tol1e-06_${BLOCKED}x${BLOCKEDCOEFF}_k1x1.0_onset24040000.0_n1_freq100.0_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1200_trainT500.0.mat does not exist"
        echo "python3 model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 1 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX $Ntrains 500 $initfile $BLOCKED $BLOCKEDCOEFF $ALTERED $ALTEREDCOEFF"
        python3 model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 1 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX $Ntrains 500 $initfile $BLOCKED $BLOCKEDCOEFF $ALTERED $ALTEREDCOEFF
      fi
    done
    BLOCKED=${GLUR},${BLOCKEDS2}
    BLOCKEDCOEFF=${GLURCOEFF},${BLOCKEDCOEFFS2[iblockcoeff]}      
    if [ -f nrn_tstop32000000_tol1e-06_${BLOCKED}x${BLOCKEDCOEFF}_k1x1.0_onset24040000.0_n1_freq100.0_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1200_trainT500.0.mat ]
    then
      echo "nrn_tstop32000000_tol1e-06_${BLOCKED}x${BLOCKEDCOEFF}_k1x1.0_onset24040000.0_n1_freq100.0_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1200_trainT500.0.mat exists"
    else
      echo "nrn_tstop32000000_tol1e-06_${BLOCKED}x${BLOCKEDCOEFF}_k1x1.0_onset24040000.0_n1_freq100.0_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1200_trainT500.0.mat does not exist"
      echo "python3 model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 1 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX $Ntrains 500 $initfile $BLOCKED $BLOCKEDCOEFF $ALTERED $ALTEREDCOEFF"
      python3 model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 1 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX $Ntrains 500 $initfile $BLOCKED $BLOCKEDCOEFF $ALTERED $ALTEREDCOEFF
    fi
  done

  #4) Simulations with combinations of protein concentrations changed
  for iblock in 10 11
  do
    BLOCKED=${GLUR},${BLOCKEDCOMBS[iblock]}
    BLOCKEDCOEFF=${GLURCOEFF},${BLOCKEDCOMBCOEFFS[iblock]}      
    if [ -f nrn_tstop32000000_tol1e-06_${GLUR},CMcombx${GLURCOEFF},${iblock}_k1x1.0_onset24040000.0_n1_freq100.0_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1200_trainT500.0.mat ]
    then
      echo "nrn_tstop32000000_tol1e-06_${GLUR},CMcombx${GLURCOEFF},${iblock}_k1x1.0_onset24040000.0_n1_freq100.0_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1200_trainT500.0.mat exists"
    else
      echo "nrn_tstop32000000_tol1e-06_${GLUR},CMcombx${GLURCOEFF},${iblock}_k1x1.0_onset24040000.0_n1_freq100.0_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1200_trainT500.0.mat does not exist"
      echo "python3 model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 1 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX $Ntrains 500 $initfile $BLOCKED $BLOCKEDCOEFF $ALTERED $ALTEREDCOEFF nrn_tstop32000000_tol1e-06_${GLUR},CMcombx${GLURCOEFF},${iblock}_k1x1.0_onset24040000.0_n1_freq100.0_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1200_trainT500.0"
            python3 model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 1 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX $Ntrains 500 $initfile $BLOCKED $BLOCKEDCOEFF $ALTERED $ALTEREDCOEFF nrn_tstop32000000_tol1e-06_${GLUR},CMcombx${GLURCOEFF},${iblock}_k1x1.0_onset24040000.0_n1_freq100.0_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1200_trainT500.0
    fi
  done
  
done
