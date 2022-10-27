function [bp] = bp_from_Chla(Chla,lambda)
% Implement the empirical relationship from Morel and Maritorena 2001 to
% calculate the particulate scattering coefficient (bp) at preselected
% output light wavelengths (lambda) using the input value of chlorophyll-a
% concentration (Chla)
%
%Reference: Morel, A., & Maritorena, S. (2001). Bio-optical properties of
%oceanic waters: A reappraisal. Journal of Geophysical Research, 106(C4),
%7163–7180. https://doi.org/10.1029/2000JC000319
%
%Inputs: Chla, lambda
%   Chla (nx1 Double): Vector containing chlorophyll-a concentrations [mg
%   m-3].
%
%   lambda (mx1 Double): Vector containing light wavelengths [nm] at which
%   output values of bp are calculated. Note that if n~=m, then either m or
%   n must be 1. Note: light wavelength is in vacuum
%
%Outputs: bp
%   bp (max(n,m)x1 Double): Particulate scattering coefficient [m-1] at
%   input light wavelengths.
%
%Created: May 13, 2022 
%Completed: August 1, 2022
%Updates: August 1, 2022
%
%Aster Taylor; modified from Sean Hart and Matthew Kehrli
%Ocean Optics Research Laboratory, Scripps Institiution of Oceanography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Particulate scattering coefficient at 660 nm from Morel and Maritorena
%(2001)
bp_660 = 0.347*Chla.^(0.766); 

%Calculation of bp at a particular input wavelength lambda
bp = bp_660.*(lambda./660).^(-1); 
end
