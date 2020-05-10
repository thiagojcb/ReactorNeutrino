# ReactorNeutrino

Scripts for reactor neutrino analysis (need ROOT to run).

## JUNO Solar Neutrino Parameters Analysis

Four main steps to generate the inputs related to signal (IBDs) for this analysis:

### 1) Expected signal
Under `JUNO_Flux_cov_matrix` directory, this step is to generate the expected number of events per energy bin (true, ie, neutrino energy, MeV). It takes two possible inputs: Hubber-Haag or DYB emission spctra. To run it, go into the respective directory, run root and do the following example commands:
```
.L GetPrediction_vDCDB_HuberHaag.cpp
MakeCovMat(10000)
```
This will generate a root file as output containing the expected, non-oscillated, spectrum and the associated error matrix, calculated with 10000 random throws.

### 2) Event tree
Under `IBD_Gen/JUNO_SNiPER_Simple` directory, this step is to generate an event tree containing the neutrino energy, the visible energy (in PE, based on SNiPER simulated response) and the event baseline. To run it, go into the respective directory, run root and do the following example commands:
```
.L IBD_gen_sPMT_SNiPER_response_simple.cpp
IBD_Gen(2020, 1e7, true)
```
This will generate a root file as output containing the event tree.
- `2020` is a date value, to set the output file name
- `1e7` is the number of events to be generated
- `true` refers to the use of the DYB spectrum. Set `false` if Huber-Haag is desirable

This macro uses the output of step 1) as input. The file name is hard-coded. So if other value instead of `10000` is used in step 1), need to update the code here accordingly.

### 3) Covariance matrix in visilbe energy (PE)
Also under `JUNO_Flux_cov_matrix` directory, this step is to generate the flux covariance matrix in visible energy (PE). To run it, go into the respective directory, run root and do the following example commands:
```
.L Convert_Eth_to_Evis_matrix.cpp
Convert_Eth_to_Evis_matrix(10000,true)
```
This will generate a root file as output containing the covariance matrix.
- `10000` is the number of throws to calculate the matrix
- `true` refers to the use of the DYB spectrum. Set `false` if Huber-Haag is desirable

This macro uses the output of step 1) and 2). The file names are hard-coded, so if `10000` of step 1) and/or `1e7` of step 2) is changed, the code need to be updated accordingly. Histogram binning is hard-coded here. If binning scheme changes, code need to be updated.

### 4) Energy Scale covariance model (DC based)
Under `JUNO_Escale_cov_matrix` directory, this step is to generate the covariance matrix associated with the energy model assumed (DC based here).
To run it, go into the respective directory, run root and do the following example commands:
```
.L MakeEscaleCovMat_DC.cpp
MakeEscaleCovMat_DC(10000, "../data/events_SPMT_1e+07evt_2020_DYB.root", true, true, false)
```
This will generate a root file as output containing the energy scale covariance matrix.
- `10000` is the number of throws to calculate the matrix
- `../data/events_SPMT_1e+07evt_2020_DYB.root` is the input file, generated at step 2) above.
- `true` refers to the use of the DYB spectrum. Set `false` if Huber-Haag is desirable. This affects only the output file name
- `true` if oscillation efects should be used (recommended). `false` to NOT consider oscillation
- `false` refering to the usual DC model to be considered. `true` if an improved model is to be considered.