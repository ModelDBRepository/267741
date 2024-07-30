:Reference :Colbert and Pan 2002

NEURON	{
	SUFFIX NaTs2_t
	USEION na READ ena WRITE ina
	RANGE gNaTs2_tbar, gNaTs2_t, ina, offm, offh, slom, sloh, tauma, taumb, tauha, tauhb
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gNaTs2_tbar = 0.00001 (S/cm2)
        offm = -32 (mV)
        offh = -60 (mV)
        slom = 6.0 (mV)
        sloh = 6.0 (mV)
        tauma = 5.49451 (ms)
        taumb = 8.06452 (ms)
        tauha = 66.6667 (ms)
        tauhb = 66.6667 (ms)
}

ASSIGNED	{
	v	(mV)
	ena	(mV)
	ina	(mA/cm2)
	gNaTs2_t	(S/cm2)
	mInf
	mTau
	mAlpha
	mBeta
	hInf
	hTau
	hAlpha
	hBeta
}

STATE	{
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gNaTs2_t = gNaTs2_tbar*m*m*m*h
	ina = gNaTs2_t*(v-ena)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
	rates()
	m = mInf
	h = hInf
}

PROCEDURE rates(){
  LOCAL qt
  qt = 2.3^((34-21)/10)
	
  UNITSOFF
    if(v == offm){
    	v = v+0.0001
    }
		mAlpha = -(offm-v)/(1-(exp((offm-v)/slom)))/tauma
		mBeta  = (offm-v)/(1-(exp(-(offm-v)/slom)))/taumb
		mTau = (1/(mAlpha + mBeta))/qt
		mInf = mAlpha/(mAlpha + mBeta)

    if(v == offh){
      v = v + 0.0001
    }

		hAlpha = (offh-v)/(1-(exp(-(offh-v)/sloh)))/tauha
		hBeta  = -(offh-v)/(1-(exp((offh-v)/sloh)))/tauhb
		hTau = (1/(hAlpha + hBeta))/qt
		hInf = hAlpha/(hAlpha + hBeta)
	UNITSON
}