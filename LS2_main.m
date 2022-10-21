function [a,anw,bb,bbp,kappa] = LS2_main(sza,lambda,Rrs,Kd,aw,bw,bp,LS2_LUT,Flag_Raman)
%Implements the LS2 inversion model to calculate a, anw, bb, and bbp from
%Rrs at specified input light wavelength
%
%Reference: Loisel, H., D. Stramski, D. Dessaily, C. Jamet, L. Li, and R.
%A. Reynolds. 2018. An inverse model for estimating the optical absorption
%and backscattering coefficients of seawater from remote-sensing
%reflectance over a broad range of oceanic and coastal marine environments.
%Journal of Geophysical Research: Oceans, 123, 2141–2171. doi:
%10.1002/2017JC013632 Note: Labeled "Steps" in this code refer to Table 1
%of this paper
%
%Required Function Inputs:
%   sza [1x1 double]: Solar zenith angle [deg]
%
%   lambda [1x1 double]: Light wavelength in vacuum [nm]; valid range is
%   302-702 nm
%
%   Rrs [1x1 double]: Spectral remote-sensing reflectance [sr^-1] at light wavelength lambda 
%
%   Kd [1x1 double]: Spectral attenuation coefficient of
%   downwelling planar irradiance, <Kd>_1 [m^-1] at lambda, averaged between the sea surface and first attenuation
%   depth  
%   
%   aw [1x1 double]: Spectral pure seawater absorption coefficient [m^-1] at
%   lambda 
%
%   bw [1x1 double]: Spectral pure seawater scattering coefficient [m^-1] at lambda 
%
%   bp [1x1 double]: Spectral particulate scattering coefficient [m^-1] at lambda
%   
%
%   LS2_LUT [1x1 struct]: Structure containing five required look-up
%   tables; can be loaded via load('LS2_LUT.mat');
%
%       LUT.muw [8x1 double]:  8 values of muw (in descending order) used
%       to construct the a and bb LUTs
%       where muw is the cosine of the angle of refraction of the solar beam just beneath the sea surface
%       
%       LUT.eta [21x1 double]: 21 values of eta (in ascending order) used 
%       to construct the a and bb LUTs
%       where eta is the ratio of the pure seawater (molecular) scattering coefficient to the total scattering coefficient of seawater
%       
%       LUT.a [21x8x4 double]: Look-up table of four polynomial
%       coefficients used to calculate the absorption coefficient from Eq.
%       9 of Loisel et al. 2018
%
%       LUT.bb  [21x8x3 double]: Look-up table of three polynomial
%       coefficients used to calculate the backscattering coefficient from
%       Eq. 8 of Loisel et al. 2018
%
%       LUT.kappa [101x7 double]: Look-up table of polynomial coefficients
%       used to calculate the dimensionless Raman scattering correction
%       factor, kappa, from the ratio bb/a. Columns represent values for
%       lambda, four polynomial coefficients (A,B,C,D), Minimum valid bb/a,
%       Maximum valid bb/a
%
%   Flag_Raman [1x1 Double]: Flag to apply or omit a Raman scattering
%   correction to Rrs If input value = 1, a Raman scattering correction is
%   applied to Rrs and output is recalculated via a single iteration. If
%   input value is not equal to 1, no Raman scattering correction is 
%   applied to Rrs and initial model output is returned
%
%Outputs:
%   a [1x1 Double]: Spectral absorption coefficient [m^-1] at lambda 
% 
%   anw [1x1 Double]: Spectral nonwater absorption coefficient [m^-1] at lambda
%   
%     
%   bb [1x1 Double]: Spectral backscattering coefficient [m^-1] at lambda
%   
%
%   bbp (1x1 Double): Spectral particulate backscattering coefficient [m^-1] at
%   lambda 
%
%   kappa [1x1 Double]: Value of the Raman scattering correction factor,
%   kappa (dim), applied to input Rrs
%
%
%Version History: 
%2018-04-04: Original implementation in C written by David Dessailly
%2020-03-23: Original Matlab version, D. Jorge 
%2022-09-01: Revised Matlab version, M. Kehrli
%DARIUSZ: WHEN WE ARE FINISHED, ADD DATE AND Final Revised Matab version: M. Kehrli, R. A. Reynolds and D. Stramski
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check function arguments and existence of LUTs
    arguments
        sza (1,1) double
        lambda (1,1) double
        Rrs (1,1) double
        Kd (1,1) double 
        aw (1,1) double
        bw (1,1) double
        bp (1,1) double
        LS2_LUT (1,1) struct
        Flag_Raman (1,1) double
    end   

%% Step 1: Calculation of muw, the cosine of the angle of refraction of the solar beam just beneath the sea surface
    nw = 1.34; %Refractive index of seawater
    muw = cosd(asind(sind(sza)/nw)); %[dim]

%% Step 2: Calculation of <Kd>_1, the average spectral attenuation coefficient of downwelling planar irradiance between the sea surface and first attenuation depth
    %In this version of the code, <Kd>_1 is assumed to be known and
    %provided as input in units of [m^-1]. In the 2018 paper, it is
    %obtained from a separate neural network algorithm

%% Step 3: Calculation of b, the total scattering coefficient in units of [m^-1]
    %In this version of the code, bp and bw are assumed to be known and
    %provided as input in units of [m^-1]. In the 2018 paper, bp is
    %estimated from chlorophyll-a concentration (Chl) where Chla is calculated from spectral remote-sensiung reflectance using the ocean color algorithm OC4v4
    b = bp + bw; %[m^-1]
%   DARIUSZ: I HAVE RESERVATIONS ABOUT STEP 3 IN A SENSE THAT THIS VERSION OF OUR CODE DOES NOT PROVIDE THE USER WITH A FULL CAPABILITY TO CALCULATE a AND bb FROM Rrs.
%   ALTHOUGH GETTING bp FROM Chl IS ONE POTENTIAL OPTION, WE DO NOT HAVE ANOTHER OPTION TO OFFER AT THIS TIME, SO I THINK WE SHOULD OFFER WITH THIS CODE A FUNCTION
%   WHICH CALCULATES bp FROM Chl, ALBEIT I AM NOT SURE IF WE SHOULD ALSO INCLUDE IN THIS OPTION THE CALCULATION OF Chl FROM REFLECTANCE SUCH AS OC4v4.
%   AS A MINIMUM, HOWEVER, I THINK WE SHOULD CONSIDER HAVING A SIMPLE FUNCTION (PROBABLY EXTERNAL TO THE MAIN CODE) TO CALCULATE bp FROM Chl. THIS WAY THE USER 
%   WILL ONLY HAVE TO WORRY ABOUT PROVIDING Chl AS ADDITIONAL INPUT TO OUR MODEL RATHER THAN ALSO HAVE A CODE THAT CONVERTS Chl TO bp.

%% Step 4: Calculation of eta, the ratio of the pure seawater (molecular) scattering coefficient to the total scattering coefficient
    eta = bw/b; %[dim]

%% Steps 5 & 7: Calculation of a and bb from Eqs. 9 and 8
    if ~isnan(eta) && ~isnan(muw)
        %find leftmost index of eta and mu values in the LUTs for
        %interpolation
        idx_eta = LS2_seek_pos(eta,LS2_LUT.eta,'eta');
        idx_muw = LS2_seek_pos(muw,LS2_LUT.muw,'muw');
        
        %If eta or mu is outside of the bounds of LUTs return nan outputs
        if isnan(idx_eta) || isnan(idx_muw)
            a = nan;
            anw = nan;
            bb = nan;
            bbp = nan;
            kappa = nan;
            return
        end

        %calculation of a from Eq. 9
        
        %a is first calculated for four combinations of LUT values
        %bracketing the lower and upper limits of eta and muw
        a00 = Kd./(LS2_LUT.a(idx_eta,idx_muw,1) + ...
            LS2_LUT.a(idx_eta,idx_muw,2).*Rrs + ...
            LS2_LUT.a(idx_eta,idx_muw,3).*Rrs.^2 + ...
            LS2_LUT.a(idx_eta,idx_muw,4).*Rrs.^3);

        a01 = Kd./(LS2_LUT.a(idx_eta,idx_muw+1,1) + ...
            LS2_LUT.a(idx_eta,idx_muw+1,2).*Rrs + ...
            LS2_LUT.a(idx_eta,idx_muw+1,3).*Rrs.^2 ...
            + LS2_LUT.a(idx_eta,idx_muw+1,4).*Rrs.^3);

        a10 = Kd./(LS2_LUT.a(idx_eta+1,idx_muw,1) + ...
            LS2_LUT.a(idx_eta+1,idx_muw,2).*Rrs + ...
            LS2_LUT.a(idx_eta+1,idx_muw,3).*Rrs.^2 + ...
            LS2_LUT.a(idx_eta+1,idx_muw,4).*Rrs.^3);

        a11 = Kd./(LS2_LUT.a(idx_eta+1,idx_muw+1,1) + ...
          LS2_LUT.a(idx_eta+1,idx_muw+1,2).*Rrs + ...
            LS2_LUT.a(idx_eta+1,idx_muw+1,3).*Rrs.^2 + ...
            LS2_LUT.a(idx_eta+1,idx_muw+1,4).*Rrs.^3);
                
        %Calculate a using 2-D linear interpolation determined from
        %bracketed values of eta and muw
        a = interp2(LS2_LUT.eta(idx_eta:idx_eta+1), ...
          LS2_LUT.muw(idx_muw:idx_muw+1), [a00 a10; a01 a11], ...
          eta, muw); %[m^-1]   

        %calculation of bb from Eq. 8
        
        %bb is first calculated for four combinations of LUT values
        %bracketing the lower and upper limits of eta and muw
        bb00 = Kd.*(LS2_LUT.bb(idx_eta,idx_muw,1).*Rrs + ...
               LS2_LUT.bb(idx_eta,idx_muw,2).*Rrs.^2 + ...
               LS2_LUT.bb(idx_eta,idx_muw,3).*Rrs.^3);

        bb01 = Kd.*(LS2_LUT.bb(idx_eta,idx_muw+1,1).*Rrs + ...
               LS2_LUT.bb(idx_eta,idx_muw+1,2).*Rrs.^2 + ...
               LS2_LUT.bb(idx_eta,idx_muw+1,3).*Rrs.^3);

        bb10 = Kd.*(LS2_LUT.bb(idx_eta+1,idx_muw,1).*Rrs + ...
               LS2_LUT.bb(idx_eta+1,idx_muw,2).*Rrs.^2 + ...
               LS2_LUT.bb(idx_eta+1,idx_muw,3).*Rrs.^3);

        bb11 = Kd.*(LS2_LUT.bb(idx_eta+1,idx_muw+1,1).*Rrs + ...
               LS2_LUT.bb(idx_eta+1,idx_muw+1,2).*Rrs.^2 + ...
               LS2_LUT.bb(idx_eta+1,idx_muw+1,3).*Rrs.^3);

        %Calculate bb using 2-D linear interpolation determined from
        %bracketed values of eta and muw
        bb = interp2(LS2_LUT.eta(idx_eta:idx_eta+1), ...
              LS2_LUT.muw(idx_muw:idx_muw+1), [bb00 bb10; bb01 bb11], ...
              eta, muw); %[m^-1]     
          
    %return NaNs if either eta or muw are NaNs
    else
        a = nan;
        anw = nan;
        bb = nan;
        bbp = nan;
        kappa = nan;
        return
    end

%% Step 9: Application of the Raman scattering correction if selected
    %If Flag_Raman is set to 1 (true), apply Raman correction to input Rrs
    %and recalculate a and bb the same as above. Otherwise no correction is
    %applied and original values are returned with kappa value of 1
    if Flag_Raman 
        %Calls subfunction LS2_calc_kappa
        kappa = LS2_calc_kappa(bb/a,lambda,LS2_LUT.kappa);

        %Apply Raman scattering correction to Rrs and recalculate a & bb
        if ~isnan(kappa)
            Rrs = Rrs.*kappa;
            
            %calculation of a from Eq. 9

            %a is first calculated for four combinations of LUT values
            %bracketing the lower and upper limits of eta and muw
            a00 = Kd./(LS2_LUT.a(idx_eta,idx_muw,1) + ...
                LS2_LUT.a(idx_eta,idx_muw,2).*Rrs + ...
                LS2_LUT.a(idx_eta,idx_muw,3).*Rrs.^2 + ...
                LS2_LUT.a(idx_eta,idx_muw,4).*Rrs.^3);

            a01 = Kd./(LS2_LUT.a(idx_eta,idx_muw+1,1) + ...
                LS2_LUT.a(idx_eta,idx_muw+1,2).*Rrs + ...
                LS2_LUT.a(idx_eta,idx_muw+1,3).*Rrs.^2 ...
                + LS2_LUT.a(idx_eta,idx_muw+1,4).*Rrs.^3);

            a10 = Kd./(LS2_LUT.a(idx_eta+1,idx_muw,1) + ...
                LS2_LUT.a(idx_eta+1,idx_muw,2).*Rrs + ...
                LS2_LUT.a(idx_eta+1,idx_muw,3).*Rrs.^2 + ...
                LS2_LUT.a(idx_eta+1,idx_muw,4).*Rrs.^3);

            a11 = Kd./(LS2_LUT.a(idx_eta+1,idx_muw+1,1) + ...
              LS2_LUT.a(idx_eta+1,idx_muw+1,2).*Rrs + ...
                LS2_LUT.a(idx_eta+1,idx_muw+1,3).*Rrs.^2 + ...
                LS2_LUT.a(idx_eta+1,idx_muw+1,4).*Rrs.^3);

            %Calculate a using 2-D linear interpolation determined from
            %bracketed values of eta and muw
            a = interp2(LS2_LUT.eta(idx_eta:idx_eta+1), ...
              LS2_LUT.muw(idx_muw:idx_muw+1), [a00 a10; a01 a11], ...
              eta, muw); %[m^-1]   

            %calculation of bb from Eq. 8

            %bb is first calculated for four combinations of LUT values
            %bracketing the lower and upper limits of eta and muw
            bb00 = Kd.*(LS2_LUT.bb(idx_eta,idx_muw,1).*Rrs + ...
                   LS2_LUT.bb(idx_eta,idx_muw,2).*Rrs.^2 + ...
                   LS2_LUT.bb(idx_eta,idx_muw,3).*Rrs.^3);

            bb01 = Kd.*(LS2_LUT.bb(idx_eta,idx_muw+1,1).*Rrs + ...
                   LS2_LUT.bb(idx_eta,idx_muw+1,2).*Rrs.^2 + ...
                   LS2_LUT.bb(idx_eta,idx_muw+1,3).*Rrs.^3);

            bb10 = Kd.*(LS2_LUT.bb(idx_eta+1,idx_muw,1).*Rrs + ...
                   LS2_LUT.bb(idx_eta+1,idx_muw,2).*Rrs.^2 + ...
                   LS2_LUT.bb(idx_eta+1,idx_muw,3).*Rrs.^3);

            bb11 = Kd.*(LS2_LUT.bb(idx_eta+1,idx_muw+1,1).*Rrs + ...
                   LS2_LUT.bb(idx_eta+1,idx_muw+1,2).*Rrs.^2 + ...
                   LS2_LUT.bb(idx_eta+1,idx_muw+1,3).*Rrs.^3);

            %Calculate bb using 2-D linear interpolation determined from
            %bracketed values of eta and muw
            bb = interp2(LS2_LUT.eta(idx_eta:idx_eta+1), ...
                  LS2_LUT.muw(idx_muw:idx_muw+1), [bb00 bb10; bb01 bb11], ...
                  eta, muw); %[m^-1]   
        end
    
    %if Flag is not 1, do nothing and return original values of a, anw, bb,
    %bbp with kappa returned as 1
    else
        kappa = 1;
    
    end

%% Steps 6 & 8: Calculation of anw and bbp
    anw = a - aw; % spectral nonwater absorption coefficient [m^-1]
    bbw = bw/2; %spectral put seawater backscattering coefficient [m^-1]
    bbp = bb - bbw; %spectral particulate backscattering coefficient [m^-1]

%% If output coefficients are negative, replace with NaN
    if a < 0 
        warning('Solution for a is negative. Output a set to nan.')
        a = nan;
    end
    if anw < 0 
        warning('Solution for anw is negative. Output anw set to nan.')
        anw = nan;
    end
    if bb < 0 
        warning('Solution for bb is negative. Output bb set to nan.')
        anw = nan;
    end
    if bbp < 0 
        warning('Solution for bbp is negative. Output bbp set to nan.')
        anw = nan;
    end

end
%end of main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Additional subfunctions that are called
function idx = LS2_seek_pos(param,LUT,type)
%MK remake of LS2 seek_pos subroutine. The subroutine finds the leftmost
%position of the input parameter in relation to its input LUT.
%
%Inputs: param, lut, type
%   param (1x1 Double): Input muw or eta value
%
%   LUT (nx1 Double): Look-up table of muw or eta values used to determine
%   coefficients in Loisel et al. 2018. If the input is associated with muw
%   the LUT must be 8x1 and sorted in descending order, and if the input is
%   associated with eta the LUT must be 21x1 and sorted in ascending order.
%
%   type (String): Characterize param input. Valid values are 'muw' or
%   'eta'. Other inputs will produce an error.
%
%Output: idx
%   idx (1x1 Double): Leftmost position/index of input param in relation to
%   its LUT.
%
%Created: September 8, 2021
%Completed: September 8, 2021
%Updates: October 15, 2022 - Added warning messages and changed logic for
%output
%
%Matthew Kehrli
%Ocean Optics Research Laboratory, Scripps Institution of Oceanography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Check function input arguments
arguments 
    param (1,1) double
    LUT (:,1) double
    type (1,1) string
end

%muw lookup
if strcmp(type,'muw')
    if length(LUT) ~= 8 || ~issorted(LUT,'descend')
        error(['Look-up table for mu_w must be a 8x1 array sorted in '...
            'descending order'])
    end
      
   if param < min(LUT)
        warning(['muw is outside the lower bound of look-up table.' ...
            ' Solutions of a and bb are output as nan.'])
       idx = nan;
   end
   if param > max(LUT)
        warning(['muw is outside the upper bound of look-up table.' ...
            ' Solutions of a and bb are output as nan.'])
       idx = nan;
   end
   for i = 1:length(LUT) - 1
       if param <= LUT(i) && param > LUT(i+1)
           idx = i;
       end
   end
   
%Eta lookup
elseif strcmp(type,'eta')
    if length(LUT) ~= 21 || ~issorted(LUT)
        error(['Look-up table for eta must be a 21x1 array sorted in '...
            'ascending order'])
    end
    if param < min(LUT)
        warning(['eta is outside the lower bound of look-up table.' ...
            ' Solutions of a and bb are output as nan.'])
        idx = nan;
    end
    if param > max(LUT)
        warning(['eta is outside the upper bound of look-up table.' ...
            ' Solutions of a and bb are output nan.'])
        idx = nan;
    end
    for i = 1:length(LUT) - 1
        if param >= LUT(i) && param < LUT(i+1)
            idx = i;
        end
    end 
end
end

function kappa = LS2_calc_kappa(bb_a,lam,rLUT)
    %The subroutine determines kappa using a linear
    %interpolation/extrapolation from the Raman scattering look-up tables.
    %
    %Inputs: bb_a, lam, r
    %   bb_a (1x1 Double): Backscattering to absorption coefficient ratio
    %   output from LS2 model.
    %
    %   lam (1x1 Double): Wavelength of bb/a ratio.
    %
    %   rLUT (101x7 Double): Look-up table for kappa.
    %
    %Output: a
    %   kappa (1x1 Double): Value of Raman scattering correction at a
    %   particular bb/a ratio and wavelength.
    %
    %Created: September 14, 2021
    %Completed: September 14, 2021
    %Updates:
    %
    %Matthew Kehrli
    %Ocean Optics Research Laboratory, Scripps Institution of Oceanography
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Check function input arguments
    arguments
        bb_a (1,1) double
        lam (1,1) double
        rLUT (101,7) double
    end
    
    %Interpolate minimum and maximum allowable bb/a values that calculate
    %kappa to the input wavelegnth.
    mins = interp1(rLUT(:,1),rLUT(:,6),lam,'linear');
    maxs = interp1(rLUT(:,1),rLUT(:,7),lam,'linear');
    
    %Check if bb/a ratio falls within min/max range
    if bb_a >= mins && bb_a <= maxs
        %Calculate all kappas for the input bb/a ratio
        kappas = rLUT(:,2).*bb_a.^3 + rLUT(:,3).*bb_a.^2 + rLUT(:,4).*bb_a ...
            + rLUT(:,5);
        %Calculate output kappa using a 1-D linear interpolation to the
        %input wavlength
        kappa = interp1(rLUT(:,1),kappas,lam,'linear');
    else
        kappa = nan;
        %Send warning message to user that no Raman correction is applied
        %to input
        warning(['No Raman Correction since bb/a value is outside of the' ...
            ' acceptable range. Kappa set to nan and no correction is' ...
            ' applied. See Raman Correction LUT.'])
    end
end
