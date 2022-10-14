# LS2
---
Loisel and Stramski 2 (LS2) is a bulk inversion model for estimating the inherent optical properites of the spectral absorption, a(λ), and backscattering b<sub>b</sub>(λ) coefficients as well as their non-water components a<sub>nw</sub>(λ) and b<sub>bp</sub>(λ) from measurments of remote sensing reflectance, R<sub>rs</sub>(λ). The model  avoids assumptions about the spectral shapes of a(λ) and b<sub>b</sub>(λ) and aims to be implemented as the first step in a mutlti-step semianalytical approach with additional partitioning models to derive optical properties of seawater constituents from both in situ and remote sensing observations of R<sub>rs</sub>(λ). The complete development and validation of the LS2 model is described in [Loisel et al. 2018](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017JC013632). In correspondance with our collaborators we have converted LS2 model source code into MATLAB file format.

This README document provides information about the files within the LS2 repository. The directories Functions, LUTs, and Test_LST on the top layer of the repository  contain a total of 14 files used to run and perform tests of the LS2 module, LS2_MK_Final.m, which are described below:

---

## LS2_MK_Final.m: 
Complete LS2 module to compute a(λ), b<sub>b</sub>(λ), a<sub>nw</sub>(λ), and b<sub>bp</sub>(λ), for *m* input samples of input R<sub>rs</sub>(λ), solar zenith angle, K<sub>d</sub>(λ), b<sub>p</sub>(λ), b<sub>w</sub>(λ), and a<sub>w</sub>(λ) at *n* input wavelength bands. See supporting documentation for further details

## Ancillary Functions
Directory of functions which trasnsform measured data into the inputs necessary for LS2. 
>**bsw_AT.m**: Computes the pure water scattering coefficient from the method of Zhang et al. 2009, for a given wavelength λ, sea temperature, and salinity. See supporting documentation for further details.\
>\
>**chl_modis_oc3m.m**: Computes the Chlorophyll-a concentration using the hybrid algorithm of Hu et al. (2012). Takes as input R<sub>rs</sub>(443), R<sub>rs</sub>(488), R<sub>rs</sub>(547), R<sub>rs</sub>(667). See supporting documentation for further details.\
>\
>**get_aw_AT.m**: Gets the pure water absorption from an LUT, only requires input wavelength. See supporting documentation for further details.\
>\
>**get_bp_AT.m**: Gets b<sub>p</sub> from Chlorophyll-a and wavelength λ. Uses the method of Loisel and Morel 1998. See supporting documentation for further details.\
>\
>**Lee_Raman_AT.m**: Computes the Raman correction using Lee et al. 2010. Takes as input a given R<sub>rs</sub>, wavelength, R<sub>rs</sub>(440), and R<sub>rs</sub>(550). See supporting documentation for further details.\
>\
>**sun_position.m**: Computes the sun position in the sky for a given latitude, longitude, and time. See supporting documentation for further details.

## Functions
Directory of functions necessary to run LS2_MK_Final.m.
>**LS2_MK_Main.m**: The main function that performs LS2 for a single input R<sub>rs</sub>(λ). The function calculates a(λ), b<sub>b</sub>(λ), a<sub>nw</sub>(λ), and b<sub>bp</sub>(λ), from input R<sub>rs</sub>(λ), solar zenith angle, K<sub>d</sub>(λ), b<sub>p</sub>(λ), b<sub>w</sub>(λ), and a<sub>w</sub>(λ) at the given R<sub>rs</sub>(λ) input wavelength. See supporting documentation for further details.\
>\
>**calc_a_MK.m**: Determine a(λ) using a linear interpolation or extrapolation from the a(λ) look-up tables associated with equation 9 in Losiel et al. 2018. See supporting documentation for further details.\
>\
>**calc_bb_MK.m**: Determine b<sub>b</sub>(λ) using a linear interpolation or extrapolation from the b<sub>b</sub>(λ) look-up tables associated with equation 8 in Losiel et al. 2018. See supporting documentation for further details.\
>\
>**calc_kappa_MK.m**: Determine $\kappa$(b<sub>b</sub>/a,λ) using a linear interpolation or extrapolation from the Raman scattering correction look-up tables. See supporting documentation for further details.\
>\
>**interp_line_MK.m**: 1-D linear interpolation or extrapolation to determine a value at a specified query point. See supporting documentation for further details.\
>\
>**restit_fit_a_MK.m**: Determine a(λ) at specific values of η and μ<sub>w</sub> within the a(λ) coefficient look-up table. See supporting documentation for further details..\
>\
>**restit_fit_bb_MK.m**: Determine b<sub>b</sub>(λ) at specific values of η and μ<sub>w</sub> within the b<sub>b</sub>(λ) coefficient look-up table. See supporting documentation for further details.\
>\
>**seek_pos_MK.m**: Finds the leftmost position of either η or μ<sub>w</sub> in relation to its look-up table. See supporting documentation for further details.

## LUTs
Directory of look-up tables (LUTs) necessary to run LS2_MK_Final.m.
>**Raman_LUT.xlsx**: LUT containing the four coefficients in the polynomial that determine the Raman correction, κ, for a given b<sub>b</sub>/a ratio as well as the minimum and maximum allowable values of b<sub>b</sub>/a at a given wavelength. These values span from 302 to 702 nm in 4 nm increments.\
>\
>**a_LUT.xlsx**: LUT containing the four coefficients that determine a(λ) in equation 9 from Loisel et al. 2018 for a given η and μ<sub>w</sub>.\
>\
>**b_LUT.xlsx**: LUT containing the three coefficients that determine b<sub>b</sub>(λ) in equation 8 from Loisel et al. 2018 for a given η and μ<sub>w</sub>.\
>\
>**aw_LUT.xlsx**: LUT containing the hyperspectral water absorption coefficient, from 180-1230 nm, at every nanometer.
>\
>**Raman_Lee_LUT.xlsx**: LUT containing coefficients for Lee et al. 2010's method for Raman correction.

## Test_LS2
Directory of functions and data used to test and evaluate the performance of LS2_MK_Final.m against the original LS2 function and measured data. The functions within this folder do not need to be in your path to run LS2_MK_Final.m.
>**Raman_RST_LS2_MK_tests.m**: Script which tests the accuracy of the version of the new version of LS2 accounting for Raman scattering against the original LS2 model and measured data from the RST dataset. The MK version of the model is tested against the original LS2 model (Figures 1 & 2) as well as against the original measured data (Figures 3 & 4) at six different wavebands associated with SeaWiFS. \
>\
>**chukchi_MODIS_LS2_AT_test.m**: Script which tests LS2 and provides an example of the use of the ancillary functions. This is run on MODIS data in the Chukchi Sea. \
>\
>**RST_LS2_AT_test.m**: Script which tests LS2 on the Remote Sensing Theory (RST) data set. These data are compiled by the OORL from in situ measurements. This scripts runs LS2 on this data set, using the ancillary functions, and compares the results to the measured values. 
>\
>**IOCCG_LS2_AT_test.m**: Script which tests LS2 on an IOCCG synthetic data set. 
>
> ### Test_Data

>Directory of data to test LS2_MK_Final.m.
>>**RST_LS2_test_MK_Raman_Nov2021.mat**: Multispectral data of 244 samples derived as a subset of the field dataset used in the 2018 paper and data outputs from the original LS2 model (run_LS2_2022.m) with a Raman scattering correction. Variables with “\_LS2” suffix were calculated by inputting Rrs and solar zenith angle into the 	original LS2 model. Note additional outputs were added to original LS2 function to retrieve a(λ), a<sub>w</sub>(λ), b<sub>b</sub>(λ), b<sub>p</sub>(λ), b<sub>w</sub>(λ), and K<sub>d</sub>(λ). The parameters a<sub>w</sub>(λ), b<sub>p</sub>(λ), b<sub>w</sub>(λ), and K<sub>d</sub>(λ) determined from the original LS2 model are used as inputs for MK’s version of LS2 (LS2_MK_NoRaman_Final.m). Variables with “\_meas” suffix are measured in situ data at six different wavebands (412, 443, 490, 510, 555, and 670 nm). \
>>\
>>**RST_LS2_data.xlsx**: An Excel document containing the RST data set. Contains location, sea surface temperature, salinity, R<sub>rs</sub>, a<sub>p</sub>, a<sub>d</sub>, a<sub>ph</sub>, a<sub>g</sub>, b<sub>b</sub>, and b<sub>bp</sub>. \
>>\
>>**IOCCG_LS2_data.xlsx**: An Excel document containing a synthetic IOCCG data set. Contains R<sub>rs</sub>, a, a<sub>w</sub>, a<sub>ph</sub>, a<sub>g</sub>, a<sub>dm</sub>, b<sub>b</sub>, b<sub>b,w</sub>, b<sub>b,ch</sub>, b<sub>b,dm</sub>, R, r<sub>rs</sub>, and K<sub>d</sub> at 0, 5, 10, 20, and 40 m of depth. Further information in the sheet itself. 
>>\
>>**AQUA_MODIS.20180801_20180831.L3m.MO.RRS.Rrs_\*.4km.nc**: R<sub>rs</sub> MODIS data, at the wavelengths replaced by the '\*'. There are 5 wavelengths, which are 443, 488, 531, 547, and 667 nm. \
>>\
>>**AQUA_MODIS.20180801_20180831.L3m.MO.SST.sst.4km.nc**: Sea surface temperature, taken from MODIS. 
>>
> ### Test_Functions
>Directory of functions to test LS2_MK_Final.m.
>>**calc_error_statistics_MK20200401.m**: Modified OORL function that calculates the statistics between measured/observed (Oi) and predicted (Pi) variables. The new function has additional output parameters compared to the original function. See supporting documentation for further details. \
>>\
>>**lsqfitm2.m**: Function written by Linhai Li to calculate the Model 2 least squares fit. See supporting documentation for further details. \
>>\
>>**parseArge.m**: Subfunction of lsqfitm2.m. See supporting documentation for further details. \
>>\

> Directory of data to test LS2_MK_Final.m.
>>**RST_LS2_test_MK_Raman_Nov2021.mat**: Multispectral data of 244 samples derived as a subset of the field dataset used in the 2018 paper and data outputs from the original LS2 model (run_LS2_2022.m) with a Raman scattering correction. Variables with “_LS2” suffix were calculated by inputting Rrs and solar zenith angle into the 	original LS2 model. Note additional outputs were added to original LS2 function to retrieve a(λ), a<sub>w</sub>(λ), b<sub>b</sub>(λ), b<sub>p</sub>(λ), b<sub>w</sub>(λ), and K<sub>d</sub>(λ). The parameters a<sub>w</sub>(λ), b<sub>p</sub>(λ), b<sub>w</sub>(λ), and K<sub>d</sub>(λ) determined from the original LS2 model are used as inputs for MK’s version of LS2 (LS2_MK_NoRaman_Final.m). Variables with “_meas” suffix are measured in situ data at six different wavebands (412, 443, 490, 510, 555, and 670 nm).
>
> ### Test_Functions
>Directory of functions to test LS2_MK_Final.m.
>>**calc_error_statistics_MK20200401.m**: Modified OORL function that calculates the statistics between measured/observed (Oi) and predicted (Pi) variables. The new function has additional output parameters compared to the original function. See supporting documentation for further details.
>>
>>**lsqfitm2.m**: Function written by Linhai Li to calculate the Model 2 least squares fit. See supporting documentation for further details.
>>
>>**parseArge.m**: Subfunction of lsqfitm2.m. See supporting documentation for further details.
>>

>>**subaxis.m**: Create axis in tiled positions. Similar to MATLAB's subplot.m. See supporting documentation for further details.
>>
> ### LS2_Output
>Directory of outputs from LS2 testing. 
>>**LS2_RST_Results.xlsx**: Output file from LS2, detailing the results of the RST dataset. \
>>\
>>**MODIS-outputs.pdf**: Results of testing on MODIS data, a large PDF file containing maps of the outputs and a spectrum averaged from 25 pixels at a specific location. \
>>\
>>**RST-outputs.pdf**: Results of testing on RST data, a large PDF file containing comparisons of the LS2 outputs versus measured values, grouped by parameter, wavelength, and sample. 
>>\
>>**IOCCG-outputs.pdf**: Results of testing on IOCCG data, a large PDF file containing comparisons of the LS2 outputs versus synthetic values, grouped by parameter, wavelength, and sample. 


---
Matthew Kehrli & Aster Taylor\
mdkehrli@ucsd.edu | mdkehrli@gmail.com\
Ocean Optics Research Laboratory, Scripps Institiution of Oceanography, University of California San Diego
