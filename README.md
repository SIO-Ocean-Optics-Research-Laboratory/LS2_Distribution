# LS2
---
Loisel and Stramski 2 (LS2) is a bulk inversion model for estimating the inherent optical properites of the spectral absorption, a(λ), and backscattering b<sub>b</sub>(λ) coefficients as well as their non-water components a<sub>nw</sub>(λ) and b<sub>bp</sub>(λ) from measurments of remote sensing reflectance, R<sub>rs</sub>(λ). The model  avoids assumptions about the spectral shapes of a(λ) and b<sub>b</sub>(λ) and aims to be implemented as the first step in a mutlti-step semianalytical approach with additional partitioning models to derive optical properties of seawater constituents from both in situ and remote sensing observations of R<sub>rs</sub>(λ). The complete development and validation of the LS2 model is described in [Loisel et al. 2018](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017JC013632). In correspondance with our collaborators we have converted LS2 model source code into MATLAB file format.

This README document provides information about the files within the LS2_Distribution repository.

---

## LS2_LUT.mat
Look-up tables (LUTs) necessary to run LS2_MK_Final.m. The structure contains five fields, each of which is a necessary to run LS2.See LS2main.m function documatation for further details about the .mat file.

## LS2_main.m
The main function that runs LS2 for a single input R<sub>rs</sub>(λ). The function calculates a(λ), b<sub>b</sub>(λ), a<sub>nw</sub>(λ), and b<sub>bp</sub>(λ), from input solar zenith angle, R<sub>rs</sub>(λ),  K<sub>d</sub>(λ),  b<sub>w</sub>(λ), a<sub>w</sub>(λ), and b<sub>p</sub>(λ), at the given R<sub>rs</sub>(λ) input wavelength. See supporting documentation for further details.

## LS2_test_run.m
Script which tests LS2 on ten samples at six different wavlength associated with SeaWiFS 

## LS2_test_run_output.xls
Spreadsheet containing the input and resulting output parameters obtained from the application of LS2 on ten samples at six different wavlengths. The file is generated by running LS2_test_run.m script.

---
Matthew Kehrli*, Aster Taylor, Rick Reynolds, & Dariusz Stramski\
*mdkehrli@ucsd.edu | mdkehrli@gmail.com\
Ocean Optics Research Laboratory, Scripps Institiution of Oceanography, University of California San Diego
