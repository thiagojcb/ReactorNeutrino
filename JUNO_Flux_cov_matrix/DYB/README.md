# Daya Bay based prediction for reactor anti-neutrino flux

The script on this directory calculates the expected reactor anti-neutrino flux based on the Daya Bay measurement.

## Inputs

The script takes in one file as input:

- `DYB_ReacCovMatrix.dat` : Contains the covariance matrix of the Daya Bay reactor antineutrino spetrum measurement ( arXiv:1607.05378 ). Hard-coded in the script can be found the binning information (eBins_DYB) and the corresponding central values (spec_DYB). This file and code values correspond to tables 12 and 13 of the arXiv pub.