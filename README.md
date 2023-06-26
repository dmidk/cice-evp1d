# cice-evp1d

# CICE evp kernels introduction

The cice-kernel project host files that have been used during the EVP refactorization and performance study conducted by Intel. The underlying refactorization idea of transforming the data representation from 2D into 1D representations were originally developed at the Danish Meteorological Institute and these transformation ideas were tested in an earlier version of the CICE code. This work builds on and extends the outcome of these tests. 

All input data required for running these benchmarks have been generated and provided by DMI in 2022. The input data is based both on the current operational setup at DMI and a setup based on the Regional Arctic Model system (RASM).
The variations of the EVP kernel implementations are all based on the upstream CICE implementation which is extracted in the file (src/bench\_v0.F90).
This repository contains the files needed for standalone teset of the 1d evp solver implemented in CICE.
