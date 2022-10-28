function [bp] = bp_from_Chla(lambda,Chla)
% Implement the empirical relationship from Morel and Maritorena 2001 to
% calculate the particulate scattering coefficient (bp) at the
% output light wavelength (lambda) using the input value of chlorophyll-a
% concentration (Chla)
%
%Reference:
%
%Morel, A., & Maritorena, S. (2001). Bio-optical properties of
%oceanic waters: A reappraisal. Journal of Geophysical Research, 106(C4),
%7163–7180. https://doi.org/10.1029/2000JC000319
%
%Inputs: Chla, lambda
%   lambda [1x1 Double]: Light wavelength [nm] at which the output value of
%   bp is calculated
%
%   Chla [1x1 Double]: Chlorophyll-a concentration [mg m-3]
%
%Outputs: bp
%   bp [1x1 Double]: Particulate scattering coefficient [m-1] at
%   input light wavelength.
%
%Created: May 13, 2022 
%Completed: August 1, 2022
%Updates: August 1, 2022
%
%Aster Taylor; modified from Sean Hart and Matthew Kehrli
%Ocean Optics Research Laboratory, Scripps Institiution of Oceanography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check function arguments and existence of LUTs
arguments
    lambda (1,1) double
    Chla (1,1) double
end   

%Particulate scattering coefficient at 660 nm from Morel and Maritorena
%(2001)
bp_660 = 0.347*Chla^(0.766); 

%Calculation of bp at a particular input wavelength lambda
bp = bp_660*(lambda./660)^(-1); 
end
