# LS2_Distribution
---

LS2 is a bulk inversion model for estimating the seawater inherent optical properties of the spectral absorption, a(λ), and backscattering, b<sub>b</sub>(λ), coefficients as well as their non-water components a<sub>nw</sub>(λ) and b<sub>bp</sub>(λ) from measurements of remote-sensing reflectance, R<sub>rs</sub>(λ). The model avoids assumptions about the spectral shapes of a(λ) and b<sub>b</sub>(λ) and can be implemented as the first step in a multi-step semianalytical approach with additional absorption partitioning models to derive absorption coefficients of seawater constituents either from in situ or remote-sensing observations of R<sub>rs</sub>(λ). The complete development and validation of the LS2 model is described in [Loisel et al. 2018] (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017JC013632). The LS2 model builds upon and extends the previous version of inverse reflectance model described in Loisel and Stramski 2000 (https://doi.org/10.1364/AO.39.003001). The LS2 model code was originally written by our collaborators from France led by Dr. Hubert Loisel. The presented version of LS2 model source code is in MATLAB file format. In this version of LS2 code we introduced a number of changes with a primary purpose to streamline the structure of the code and facilitate its application by users.

This README document provides information about the files within the LS2_Distribution repository.

---

## LS2_main.m
The main function that runs LS2 for a single input R<sub>rs</sub>(λ). The function calculates a(λ), b<sub>b</sub>(λ), a<sub>nw</sub>(λ), and b<sub>bp</sub>(λ), from input solar zenith angle, R<sub>rs</sub>(λ),  K<sub>d</sub>(λ),  b<sub>w</sub>(λ), a<sub>w</sub>(λ), and b<sub>p</sub>(λ), at the given R<sub>rs</sub>(λ) input wavelength. See supporting documentation for further details.## LS2_main.m
The main function that runs LS2 for a single input R<sub>rs</sub>(λ). The function calculates a(λ), b<sub>b</sub>(λ), a<sub>nw</sub>(λ), and b<sub>bp</sub>(λ), from input solar zenith angle, R<sub>rs</sub>(λ),  K<sub>d</sub>(λ),  b<sub>w</sub>(λ), a<sub>w</sub>(λ), and b<sub>p</sub>(λ), at the given R<sub>rs</sub>(λ) input wavelength. See supporting documentation for further details.

## LS2_LUT.mat
Look-up tables (LUTs) necessary to run LS2_main.m. The structure contains five fields, each of which is necessary to run LS2. See LS2_main.m function documentation for further details about the .mat file.

## LS2_test_run.m
Script which tests LS2 on ten sample inputs at six different wavelengths corresponding to center wavelengths of spectral bands available on ocean color sensor SeaWiFS 

## LS2_test_run.xls
Excel spreadsheet containing the input and resulting output parameters obtained from the application of LS2 on ten samples at six different wavlengths. The file is the original output file generated by LS2_test_run.m.

## bp_from_Chla.m
Function that calculates the particulate scattering coefficient b<sub>p</sub>(λ) at a user defined wavelength from the emperical relationship described in Morel and Maritorena 2001. This algorithm is used to produce input of b<sub>p</sub>(λ) to the LS2 function and is provided for the convenience of the user.

---
Contributors: Matthew Kehrli, Aster Taylor, Rick A. Reynolds, and Dariusz Stramski\
Contacts: Matthew Kehrli (mdkehrli@ucsd.edu | mdkehrli@gmail.com), Rick Reynolds (rreynolds@ucsd.edu), Dariusz Stramski (dstramski@ucsd.edu)\
Ocean Optics Research Laboratory, Scripps Institution of Oceanography, University of California San Diego
