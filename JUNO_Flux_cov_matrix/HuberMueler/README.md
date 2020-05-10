# Huber-Mueler-Haag based prediction for reactor anti-neutrino flux

The script on this directory calculates the expected reactor anti-neutrino flux based on the Huber-Mueler-Haag prediction.

## Inputs

The script takes in 3 types of inputs:

- `gm_xsec.dat` : It the IBD cross-section shape values. The first column is for the neutrino eneryg (in MeV) and the second column is the cross-section value. The cross-section need to be normalized by a normalization (K) factor, that currently is hard-coded.

- `*_v1.txt` : these files contain each fissile isotope anti-neutrino emission spectra. The first column is for the energy bin lower edge (MeV), the second column for the energy bin upper edge (MeV), the third column is the number of antineutrino emitted per fission per MeV, the fourth column is the total error (in %) of this estimation. The last two columns if for the correlated and uncorrleated uncertainty (in %) between the isotopes.

- `FluxMatrixFromDCDB.root` : contains the covariance and correlation matrices of the neutrino emission prediction of the four isotopes. The order is: 235U - 238U - 239Pu - 241Pu