# Simple IBD event generator

The script on this directory generates IBD events (positrons, in MeV or PE units) based on the flux prediciton of ther scripts and the JUNO detector response simulation.

## Inputs

The script takes in one file as input:

- `PE2MeV_vs_MeV_vs_R3.txt` : This file contains the PE to MeV values given by the SNiPER code of the JUNO detector simulation. To take into account non-linearities (light and geometry) the PE to MeV was retrieved for 20 values of mean energy and 20 values of mean volume. The columns in this file represents:

Energy Counter, Volume Counter, PE2MeV_cv, PE2MeV_sigma, E_resol_%, gaus prob, nEvents on original histo

Lowest energy value: 0.5 MeV
Highest energy value: 10.5 MeV

Lowest volume value: 0 m^3
Highest volume value: 5400 m^3