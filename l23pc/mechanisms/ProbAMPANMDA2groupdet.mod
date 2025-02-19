TITLE AMPA and NMDA receptor with presynaptic short-term plasticity 


COMMENT
AMPA and NMDA receptor conductance using a dual-exponential profile
presynaptic short-term plasticity based on Fuhrmann et al. 2002
Implemented by Srikanth Ramaswamy, Blue Brain Project, July 2009
Etay: changed weight to be equal for NMDA and AMPA, gmax accessible in Neuron
Tuomo: Added the group property - each synapse may correspond to a group of synapses.
       The activated synapses are determined in advance by setVec2().
ENDCOMMENT


NEURON {

        POINT_PROCESS ProbAMPANMDA2groupdet
        RANGE tau_r_AMPA, tau_d_AMPA, tau_r_NMDA, tau_d_NMDA, Nsyns, Nevents, eventCounter
        RANGE Use, u, Dep, Fac, u0, weight_factor_NMDA
        RANGE i, i_AMPA, i_NMDA, g_AMPA, g_NMDA, e, gmax
        NONSPECIFIC_CURRENT i, i_AMPA,i_NMDA
	POINTER rng
}

PARAMETER {

        tau_r_AMPA = 0.2   (ms)  : dual-exponential conductance profile
        tau_d_AMPA = 1.7    (ms)  : IMPORTANT: tau_r < tau_d
	tau_r_NMDA = 0.29   (ms) : dual-exponential conductance profile
        tau_d_NMDA = 43     (ms) : IMPORTANT: tau_r < tau_d
        Use = 1.0   (1)   : Utilization of synaptic efficacy (just initial values! Use, Dep and Fac are overwritten by BlueBuilder assigned values) 
        Dep = 100   (ms)  : relaxation time constant from depression
        Fac = 10   (ms)  :  relaxation time constant from facilitation
        e = 0     (mV)  : AMPA and NMDA reversal potential
	mg = 1   (mM)  : initial concentration of mg2+
        mggate
    	gmax = .001 (uS) : weight conversion factor (from nS to uS)
    	u0 = 0 :initial value of u, which is the running value of Use
	Nsyns = 10 : How many synapses are there actually (the size of "space" divided by three)
	Nevents = 0 : How many events will there be (the size of "space2")
        weight_factor_NMDA = 1
}

COMMENT
The Verbatim block is needed to generate random nos. from a uniform distribution between 0 and 1 
for comparison with Pr to decide whether to activate the synapse or not
ENDCOMMENT
   
VERBATIM

#include<stdlib.h>
#include<stdio.h>
#include<math.h>

double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);

extern int ifarg(int iarg);
extern int vector_capacity(void* vv);
extern void* vector_arg(int iarg);
ENDVERBATIM
  

ASSIGNED {

        v (mV)
        i (nA)
	i_AMPA (nA)
	i_NMDA (nA)
        g_AMPA (uS)
	g_NMDA (uS)
        factor_AMPA
	factor_NMDA
	rng
 	space        : A pointer to the vector containing the synapse times. Note that the underlying vector should not be touched after initialization by setVec().
 	space2       : A pointer to the vector containing the event IDs. Note that the underlying vector should not be touched after initialization by setVec2().
        eventCounter : An index for space2 (counts the passed events)
}

STATE {

        A_AMPA       : AMPA state variable to construct the dual-exponential profile - decays with conductance tau_r_AMPA
        B_AMPA       : AMPA state variable to construct the dual-exponential profile - decays with conductance tau_d_AMPA
	A_NMDA       : NMDA state variable to construct the dual-exponential profile - decays with conductance tau_r_NMDA
        B_NMDA       : NMDA state variable to construct the dual-exponential profile - decays with conductance tau_d_NMDA
}

INITIAL{

        LOCAL tp_AMPA, tp_NMDA
        
	A_AMPA = 0
        B_AMPA = 0
	
	A_NMDA = 0
	B_NMDA = 0
        
	tp_AMPA = (tau_r_AMPA*tau_d_AMPA)/(tau_d_AMPA-tau_r_AMPA)*log(tau_d_AMPA/tau_r_AMPA) :time to peak of the conductance
	tp_NMDA = (tau_r_NMDA*tau_d_NMDA)/(tau_d_NMDA-tau_r_NMDA)*log(tau_d_NMDA/tau_r_NMDA) :time to peak of the conductance
        
	factor_AMPA = -exp(-tp_AMPA/tau_r_AMPA)+exp(-tp_AMPA/tau_d_AMPA) :AMPA Normalization factor - so that when t = tp_AMPA, gsyn = gpeak
        factor_AMPA = 1/factor_AMPA
	
	factor_NMDA = -exp(-tp_NMDA/tau_r_NMDA)+exp(-tp_NMDA/tau_d_NMDA) :NMDA Normalization factor - so that when t = tp_NMDA, gsyn = gpeak
        factor_NMDA = 1/factor_NMDA
	eventCounter = 0

}

BREAKPOINT {

        SOLVE state METHOD cnexp
        mggate = 1 / (1 + (mg/4.1 (mM))*exp(0.063 (/mV)*(-v))) :mggate kinetics - Spruston et al. 1995
        g_AMPA = gmax*(B_AMPA-A_AMPA) :compute time varying conductance as the difference of state variables B_AMPA and A_AMPA
	g_NMDA = gmax*(B_NMDA-A_NMDA) * mggate :compute time varying conductance as the difference of state variables B_NMDA and A_NMDA and mggate kinetics
        i_AMPA = g_AMPA*(v-e) :compute the AMPA driving force based on the time varying conductance, membrane potential, and AMPA reversal
	i_NMDA = g_NMDA*(v-e) :compute the NMDA driving force based on the time varying conductance, membrane potential, and NMDA reversal
	i = i_AMPA + i_NMDA
}

DERIVATIVE state{

        A_AMPA' = -A_AMPA/tau_r_AMPA
        B_AMPA' = -B_AMPA/tau_d_AMPA
	A_NMDA' = -A_NMDA/tau_r_NMDA
        B_NMDA' = -B_NMDA/tau_d_NMDA
}


NET_RECEIVE (weight, Pv, Pr, u, myInd, tsyn (ms), Pv_tmp){
	
	:printf("NMDA weight = %g\n", weight_NMDA)

        INITIAL{
                Pv=1
                u=u0
		eventCounter=0
            }

        :Randomize which of the synapses is activated. Note that an additional random number is generated by rand() - this may interfere with the random number order in parallel simulations.
        VERBATIM
          void** vv = (void**)(&space);
          void** vv2 = (void**)(&space2);
          double *x;
          int nx = vector_instance_px(*vv, &x);
          double *x2;
          int nx2 = vector_instance_px(*vv2, &x2);
          int myInd = 0;
          if (eventCounter < nx2) {
            myInd = x2[(int)eventCounter];
            //printf("eventCounter < nx2. t = %lf, eventCounter = %lf, nx2 = %i, myInd = %i\n",t, eventCounter, nx2, myInd);
          }
          else printf("eventCounter >= nx2! t = %lf, eventCounter = %lf, nx2 = %i\n",t, eventCounter, nx2);
	  eventCounter++;
          _args[4] = myInd;
          _args[5] = x[myInd];                //tsyn
          _args[1] = x[myInd+(int)Nsyns];     //Pv
          _args[3] = x[myInd+2*((int)Nsyns)]; //u
        ENDVERBATIM
	::printf("NET_RECEIVE_beg: Pv = %g, Pr = %g, u = %g, myInd = %g, tsyn = %g, t = %g\n", Pv, Pr, u, myInd, tsyn, t)
	:printf("NET_RECEIVE_beg:  myInd = %g/%g, Pv = %g, u = %g, tsyn = %g, t = %g. ", myInd, Nsyns, Pv, u, tsyn, t)

        : calc u at event-
        if (Fac > 0) {
              u = u*exp(-(t - tsyn)/Fac) :update facilitation variable if Fac>0 Eq. 2 in Fuhrmann et al.
        } else {
              u = Use  
        } 
        if(Fac > 0){
              u = u + Use*(1-u) :update facilitation variable if Fac>0 Eq. 2 in Fuhrmann et al.
        }    

        
        Pv_tmp  = 1 - (1-Pv) * exp(-(t-tsyn)/Dep) :Probability Pv for a vesicle to be available for release, analogous to the pool of synaptic
                                                  :resources available for release in the deterministic model. Eq. 3 in Fuhrmann et al.
        Pr  = u * Pv_tmp                          :Pr is calculated as Pv * u (running value of Use)
        Pv_tmp  = Pv_tmp - u * Pv_tmp             :update Pv as per Eq. 3 in Fuhrmann et al.
        :printf("Pv = %g\n", Pv)
        :printf("Pr = %g\n", Pr)

	if (erand() < Pr){
            tsyn = t
	    Pv = Pv_tmp
            A_AMPA = A_AMPA + weight*factor_AMPA
            B_AMPA = B_AMPA + weight*factor_AMPA
            A_NMDA = A_NMDA + weight*weight_factor_NMDA*factor_NMDA
            B_NMDA = B_NMDA + weight*weight_factor_NMDA*factor_NMDA
	    :::printf("Exc Released! Pr = %g, Pv = %g, u = %g, myInd = %g, tsyn = %g\n", Pr, Pv, u, myInd, tsyn)
	    ::printf("Exc Released!\n")
            :printf ( "R! Pr = %g\n" , Pr )
          } else {
	    :::printf("Exc Not released! Pr = %g, Pv = %g, u = %g, myInd = %g, tsyn = %g\n", Pr, Pv, u, myInd, tsyn)
            ::printf("Exc Not released! value = %g, Pr = %g\n", erand(), Pr )
            :printf ( "NR! Pr = %g\n" , Pr )
          }
	:printf("NET_RECEIVE_end: Pv = %g, Pr = %g, u = %g, myInd = %g, tsyn = %g, t = %g\n", Pv, Pr, u, myInd, tsyn, t)

        VERBATIM
          x[myInd] = _args[5];
          x[myInd+(int)Nsyns] = _args[1];
	  x[myInd+2*((int)Nsyns)] = _args[3];
        ENDVERBATIM
}

PROCEDURE setRNG() {
VERBATIM
    {
        /**
         * This function takes a NEURON Random object declared in hoc and makes it usable by this mod file.
         * Note that this method is taken from Brett paper as used by netstim.hoc and netstim.mod
         * which points out that the Random must be in negexp(1) mode
         */
        void** pv = (void**)(&_p_rng);
        if( ifarg(1)) {
            *pv = nrn_random_arg(1);
        } else {
            *pv = (void*)0;
        }
    }
ENDVERBATIM
}

FUNCTION erand() {
VERBATIM
	    //FILE *fi;
        double value;
        if (_p_rng) {
                /*
                :Supports separate independent but reproducible streams for
                : each instance. However, the corresponding hoc Random
                : distribution MUST be set to Random.negexp(1)
                */
                value = nrn_random_pick(_p_rng);
		        //fi = fopen("RandomStreamMCellRan4.txt", "w");
                //fprintf(fi,"random stream for this simulation = %lf\n",value);
                //printf("random stream for this simulation = %lf\n",value);
                return value;
        }else{
ENDVERBATIM
                : the old standby. Cannot use if reproducible parallel sim
                : independent of nhost or which host this instance is on
                : is desired, since each instance on this cpu draws from
                : the same stream
                erand = exprand(1)
VERBATIM
        }
ENDVERBATIM
        :erand = value :This line must have been a mistake in Hay et al.'s code, it would basically set the return value to a non-initialized double value.
                       :The reason it sometimes works could be that the memory allocated for the non-initialized happened to contain the random value
                       :previously generated. However, here we commented this line out.
}

PROCEDURE setVec() {    : Sets the times of firing of each synapse. This should be done only once for each ProbAMPANMDA2group,
                        : before the running of the simulation, and the underlying vector should be untouched after that.
  VERBATIM
  void** vv;
  vv = (void**)(&space);
  *vv = (void*)0;
  if (ifarg(1)) {
    *vv = vector_arg(1);
    Nsyns = vector_capacity(*vv)/3;
  }
  ENDVERBATIM
}
PROCEDURE setVec2() {    : Sets the IDs of the synapses to fire
  VERBATIM
  void** vv;
  vv = (void**)(&space2);
  *vv = (void*)0;
  if (ifarg(1)) {
    *vv = vector_arg(1);
    Nevents = vector_capacity(*vv);
  }
  ENDVERBATIM
}

PROCEDURE printVec() { : Prints the previous times of firing of each synapse.
VERBATIM
    void** vv = (void**)(&space);
    double *x;
    int nx = vector_instance_px(*vv, &x);
    int i1;
    for (i1=0; i1<Nsyns;i1++) {
      printf("tsyns[%i] = %g, Pv[%i] = %g, u[%i] = %g\n", i1, x[i1], i1, x[i1+(nx/3)], i1, x[i1+2*(nx/3)]);
    }
ENDVERBATIM
}
PROCEDURE printVec2() { : Prints the previous times of firing of each synapse.
VERBATIM
    void** vv = (void**)(&space2);
    double *x;
    int nx = vector_instance_px(*vv, &x);
    int i1;
    for (i1=0; i1<Nevents;i1++) {
      printf("%g ", x[i1]);
      if (i1%100==99)
        printf("\n");
    }
ENDVERBATIM
}

