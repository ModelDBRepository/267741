Python and NEURON scripts for running plasticity simulations of cortical cells as described in the article Genetic mechanisms for impaired synaptic plasticity in schizophrenia revealed by computational modelling (Maki-Marttunen et al. 2024).
The entry contains three folders:

### syn
The main folder for plasticity simulations.

### l23pc
The folder for simulations of layer II/III pyramidal cells. Contains the models
from Markram et al. 2015 "Reconstruction and simulation of neocortical microcircuitry.
*Cell*, 163(2), 456-492.", namely, models `cADpyr229_L23_PC_5ecbf9b163`,
    `cADpyr229_L23_PC_8ef1aa6602`, `cADpyr229_L23_PC_863902f300`, `cADpyr229_L23_PC_c292d67a2e`,
    and `cADpyr229_L23_PC_c2e79db05a`.

### genetic
The folder for genetic analyses.

```
Tuomo Maki-Marttunen, 2020-2023
CC BY 4.0
```

The codes have been tested running Python versions 3.9.1 - 3.10.6 and NEURON versions 7.7 - 7.8.2.