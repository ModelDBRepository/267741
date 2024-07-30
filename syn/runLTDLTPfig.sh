FREQS=(0.001 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0 21.0 22.0 23.0 24.0 25.0 26.0 27.0)
NSTIMS=(0 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500 2600 2700)

T=27000000
ONSET=24040000
CAFLUX=100.0
LFLUX=5.0
GLUFLUX=10.0
ACHFLUX=10.0
initfile=None
BLOCKEDS=('Ca' 'Calbin,CalbinC' 'Calbin,CalbinC' 'DAGK' 'DAGK' 'Gi' 'Gi' 'MGluR' 'MGluR' 'Gqabg' 'Gqabg' 'PKC' 'PKC' 'PLA2' 'PLA2' 'CK' 'CK' 'Ng' 'Ng' 'NCX' 'NCX' 'PDE4' 'PDE4' 'PKA' 'PKA' 'PP1' 'PP1')
BLOCKEDCOEFFS=(1.0 0.8,0.8 1.2,1.2 0.8 1.2 0.8 1.2 0.8 1.2 0.8 1.2 0.8 1.2 0.8 1.2 0.8 1.2 0.8 1.2 0.8 1.2 0.8 1.2 0.8 1.2 0.8 1.2)

GLUR=GluR1,GluR1_memb,GluR2,GluR2_memb
GLURCOEFF=0.5,0.5,1.5,1.5
for iblock in `seq 0 26`
do
  BLOCKED=${GLUR},${BLOCKEDS[iblock]}
  BLOCKEDCOEFF=${GLURCOEFF},${BLOCKEDCOEFFS[iblock]}

  for ifreq in `seq 0 27`
  do	       
    FREQ=${FREQS[ifreq]}
    NSTIM=${NSTIMS[ifreq]}
    if [ -f nrn_tstop27000000_tol1e-06_${BLOCKED}x${BLOCKEDCOEFF}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1_trainT1000.0.mat ]
    then
      echo "nrn_tstop27000000_tol1e-06_${BLOCKED}x${BLOCKEDCOEFF}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1_trainT1000.0.mat exists"
    else
      echo "nrn_tstop27000000_tol1e-06_${BLOCKED}x${BLOCKEDCOEFF}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1_trainT1000.0.mat does not exist"
      echo "python3 model_nrn_altered_noU_extfilename_smallconcs.py ${T} 1e-6 $ONSET $NSTIM $FREQ 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 1000 $initfile $BLOCKED $BLOCKEDCOEFF"
      python3 model_nrn_altered_noU_extfilename_smallconcs.py ${T} 1e-6 $ONSET $NSTIM $FREQ 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 1000 $initfile $BLOCKED $BLOCKEDCOEFF
    fi
  done
done


BLOCKEDS=('Calbin,CalbinC,DAGK,Gi,MGluR,NCX,Ng,PDE4,PKC,PP1' 'Calbin,CalbinC,DAGK,Gi,MGluR,NCX,Ng,PDE4,PKC,PP1' 'CK,DAGK,Gqabg,PKA' 'DAGK,PKA,PKC,PLA2,PP1' 'CK,Calbin,CalbinC,DAGK,Gi,Gqabg,MGluR,NCX,Ng,PDE4,PKA,PKC,PP1' 'Calbin,CalbinC,DAGK,Gi,MGluR,NCX,Ng,PDE4,PKA,PKC,PLA2,PP1' 'Calbin,CalbinC,DAGK,Gi,Ng,PDE4,PKC' 'DAGK,Gi,NCX,Ng,PDE4,PP1' 'CK,DAGK,Gqabg,PKA' 'DAGK,PKA,PKC,PLA2,PP1' 'CK,Calbin,CalbinC,DAGK,Gi,Gqabg,Ng,PDE4,PKA,PKC' 'DAGK,Gi,NCX,Ng,PDE4,PKA,PKC,PLA2,PP1' 'DAGK,MGluR,NCX,PDE4,PP1' 'MGluR' 'CK,DAGK,Gqabg,PKA' 'DAGK,PKA,PKC,PLA2,PP1' 'CK,DAGK,Gqabg,MGluR,NCX,PDE4,PKA,PP1' 'DAGK,MGluR,PKA,PKC,PLA2,PP1' 'DAGK,PDE4' 'CK,DAGK,Gqabg,PKA' 'DAGK,PKA,PKC,PLA2,PP1' 'CK,DAGK,Gqabg,PDE4,PKA' 'DAGK,PKA,PKC,PLA2,PP1')
BLOCKEDCOEFFS=(0.967,0.967,0.974,0.979,1.062,1.060,0.973,1.054,0.994,1.066 1.025,1.025,0.988,0.968,1.093,1.021,0.941,0.967,0.984,0.986 0.905,1.140,1.109,1.079 1.158,1.055,1.056,0.864,0.897 0.905,0.967,0.967,1.029,0.979,1.109,1.062,1.060,0.973,1.054,1.079,0.994,1.066 1.025,1.025,1.045,0.968,1.093,1.021,0.941,0.967,1.055,1.020,0.864,0.942 0.985,0.985,0.924,0.971,0.962,1.012,0.996 0.992,0.984,0.999,0.944,0.977,0.992 0.959,1.072,1.063,1.048 1.134,1.039,1.037,0.933,0.940 0.959,0.985,0.985,0.998,0.971,1.063,0.962,1.012,1.048,0.996 1.039,0.984,0.999,0.944,0.977,1.039,1.037,0.933,0.966 0.918,1.062,1.090,1.054,1.066 1.093 0.905,1.140,1.109,1.079 1.158,1.055,1.056,0.864,0.897 0.905,1.029,1.109,1.062,1.090,1.054,1.079,1.066 1.158,1.093,1.055,1.056,0.864,0.897 0.924,1.012 0.959,1.072,1.063,1.048 1.134,1.039,1.037,0.933,0.940 0.959,0.998,1.063,1.012,1.048 1.134,1.039,1.037,0.933,0.940)

for icomb in 10 11
do
  BLOCKED=${GLUR},${BLOCKEDS[icomb]}
  BLOCKEDCOEFF=${GLURCOEFF},${BLOCKEDCOEFFS[icomb]}
  for ifreq in `seq 0 27`
  do	       
    FREQ=${FREQS[ifreq]}
    NSTIM=${NSTIMS[ifreq]}
    if [ -f nrn_tstop27000000_tol1e-06_${GLUR},CMcombx${GLURCOEFF},${icomb}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1_trainT1000.0.mat ]
    then
      echo "nrn_tstop27000000_tol1e-06_${GLUR},CMcombx${GLURCOEFF},${icomb}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1_trainT1000.0.mat exists"
    else
      echo "nrn_tstop27000000_tol1e-06_${GLUR},CMcombx${GLURCOEFF},${icomb}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1_trainT1000.0.mat does not exist"
      echo "python3 model_nrn_altered_noU_extfilename_smallconcs.py ${TSHORT} 1e-6 $ONSET $NSTIM $FREQ 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 1000 $initfile $BLOCKED $BLOCKEDCOEFF 1 1.0 nrn_tstop27000000_tol1e-06_${GLUR},CMcombx${GLURCOEFF},${icomb}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1_trainT1000.0"
      python3 model_nrn_altered_noU_extfilename_smallconcs.py ${TSHORT} 1e-6 $ONSET $NSTIM $FREQ 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 1000 $initfile $BLOCKED $BLOCKEDCOEFF 1 1.0 nrn_tstop27000000_tol1e-06_${GLUR},CMcombx${GLURCOEFF},${icomb}_onset24040000.0_n${NSTIM}_freq${FREQ}_dur3.0_flux${CAFLUX}_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1_trainT1000.0
    fi
  done
done
