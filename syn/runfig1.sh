
#1) Run first one simulation of each of the three types to get baseline values. This takes 3-4 min each.

GLUR=GluR1,GluR1_memb,GluR2,GluR2_memb
GLURCOEFFS=(0.5,0.5,1.5,1.5 2.0,2.0,0.0,0.0 2.0,2.0,0.0,0.0)
BLOCKEDS=('Ca' 'PP1' 'R,PP1')
BLOCKEDCOEFFS=(1.0 0.5 0.0,0.0)
LFLUXES=(5.0 0.2 5.0)
GLUFLUXES=(10.0 10.0 10.0)
CAFLUXES=(150.0 150.0 150.0)

FREQS=(0.2 0.5 1.0 2.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 20.0 22.0 25.0 30.0 40.0 50.0 60.0 75.0 100.0 200.0 300.0 500.0)
NSTIMS=(20 50 100 200 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 2000 2200 2500 3000 4000 5000 6000 7500 10000 20000 30000 50000)


TSHORT=27000000
ONSET=24040000

GLUFLUX=10.0
ACHFLUX=10.0
initfile=None

ALTEREDS=(     1   182,185,188,203,206,209,246,249,252,267,270,273 373,374 416,417,418 )
ALTEREDCOEFFS=(1.0 0,0,0,0,0,0,0,0,0,0,0,0                         0,0     0,0,0       )

ALTEREDBLOCKEDS=(     Ca  Ca  PKAc                  PKCtCa)
ALTEREDBLOCKEDCOEFFS=(1.0 1.0 3.5544295672723423e12 20.712191443231298e12) #Baseline values from python3 printbaselinePKA_PKC.py nrn_tstop27000000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,Cax0.5,0.5,1.5,1.5,1.0_k1x1.0_onset24040000.0_n100_freq1.0_dur3.0_flux*

for iblock in 0 1 2
do
 for myiFREQ in 0
 do
  BLOCKED=${GLUR},${BLOCKEDS[iblock]}
  BLOCKEDCOEFF=${GLURCOEFFS[iblock]},${BLOCKEDCOEFFS[iblock]}
  LFLUX=${LFLUXES[iblock]}
  CAFLUX=${CAFLUXES[iblock]}
  GLUFLUX=${GLUFLUXES[iblock]}
  FREQ=${FREQS[myiFREQ]}
  NSTIM=${NSTIMS[myiFREQ]}

  ialtered=0
  MYBLOCKED=${BLOCKED},${ALTEREDBLOCKEDS[ialtered]}
  MYBLOCKEDCOEFF=${BLOCKEDCOEFF},${ALTEREDBLOCKEDCOEFFS[ialtered]}
  ALTERED=${ALTEREDS[ialtered]}
  ALTEREDCOEFF=${ALTEREDCOEFFS[ialtered]}  
  echo "python3 model_nrn_altered_noU_extfilename_smallconcs.py ${TSHORT} 1e-6 $ONSET $NSTIM $FREQ 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 1000 $initfile $MYBLOCKED $MYBLOCKEDCOEFF $ALTERED $ALTEREDCOEFF nrn_${NSTIM}_${FREQ}_${CAFLUX}_${LFLUX}_${GLUFLUX}_${ACHFLUX}_${MYBLOCKED}x${MYBLOCKEDCOEFF}_k${ALTERED}x${ALTEREDCOEFF}"
  python3 model_nrn_altered_noU_extfilename_smallconcs.py ${TSHORT} 1e-6 $ONSET $NSTIM $FREQ 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 1000 $initfile $MYBLOCKED $MYBLOCKEDCOEFF $ALTERED $ALTEREDCOEFF nrn_${NSTIM}_${FREQ}_${CAFLUX}_${LFLUX}_${GLUFLUX}_${ACHFLUX}_${MYBLOCKED}x${MYBLOCKEDCOEFF}_k${ALTERED}x${ALTEREDCOEFF}
 done
done

#2) Run iterations to find the correct reaction rate coefficients for spontaneous PKA phosphorylation:
python3 doiters_fig1.py nrn_20_0.2_150.0_5.0_10.0_10.0_GluR1,GluR1_memb,GluR2,GluR2_memb,Ca,Cax0.5,0.5,1.5,1.5,1.0,1.0_k1x1.0.mat
python3 doiters_fig1.py nrn_20_0.2_150.0_0.2_10.0_10.0_GluR1,GluR1_memb,GluR2,GluR2_memb,PP1,Cax2.0,2.0,0.0,0.0,0.5,1.0_k1x1.0.mat
python3 doiters_fig1.py nrn_20_0.2_150.0_5.0_10.0_10.0_GluR1,GluR1_memb,GluR2,GluR2_memb,R,PP1,Cax2.0,2.0,0.0,0.0,0.0,0.0,1.0_k1x1.0.mat

#3) Run the simulations with varying frequencies for control case (no spontaneous PKA or PKC phosphorylation)
python3 runvaryfreqnonaltereds_nonlog_fig1.py 0
python3 runvaryfreqnonaltereds_nonlog_fig1.py 1
python3 runvaryfreqnonaltereds_nonlog_fig1.py 2

#4) Run the simulations with varying frequencies for icase 0 and 3 (spontaneous PKA or PKC phosphorylation)
python3 runvaryfreqaltereds2_nonlog_fig1.py 0
python3 runvaryfreqaltereds2_nonlog_fig1.py 1
python3 runvaryfreqaltereds2_nonlog_fig1.py 2


