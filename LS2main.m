function [a,anw,bb,bbp,kappa] = LS2main(sza,lambda,Rrs,Kd,aw,bw,bp,LS2_LUT,Flag_Raman)
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
%   Rrs [1x1 double]: Spectral remote sensing reflectance at lambda [sr^-1]
%
%   Kd [1x1 double]: Average spectral attenuation coefficient between
%   surface and first attenuation depth for downwelling planar irradiance,
%   <Kd>_1, at lambda  [m^-1]
%   
%   aw [1x1 double]: Pure seawater spectral absorption coefficient at
%   lambda [m^-1]
%
%   bw [1x1 double]: Pure seawater scattering coefficient at lambda [m^-1]
%
%   bp [1x1 double]: Spectral particulate scattering coefficient at lambda
%   [m^-1]
%
%   LS2_LUT [1x1 struct]: Structure containing five required look-up
%   tables; can be loaded via load('LS2_LUT.mat');
%
%       LUT.muw [8x1 double]:  8 values of muw (in descending order) used
%       to construct the a and bb LUTs 
%       
%       LUT.eta [21x1 double]: 21 values of eta (in ascending order) used 
%       to construct the a and bb LUTs
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
%   applied to Rrs and output is recalculated via a single iteration If
%   input value ne 1, no Raman scattering correction is applied to Rrs and
%   initial model output is returned
%
%Outputs:
%   a [1x1 Double]: Spectral absorption coefficient at lambda [m^-1]
% 
%   anw [1x1 Double]: Spectral nonwater absorption coefficient at lambda
%   [m^-1]
%     
%   bb [1x1 Double]: Spectral backscattering coefficient at lambda
%   [m^-1]
%
%   bbp (1x1 Double): Spectral particulate backscattering coefficient at
%   lambda [m^-1]
%
%   kappa [1x1 Double]: Value of the Raman scattering correction factor,
%   kappa, applied to input Rrs (dim)
%
%
%Version History: 
%2018-04-04: Original implementation in C written by David Dessailly
%2020-03-23: Matlab version, D. Schaeffer 
%2022-09-01: Revised Matlab version, M. Kehrli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check function arguments and existence of LUTs
    arguments  %RR - NOT SURE WHAT THIS DOES IF THERE IS AN ERROR
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

%% Step 2: Calculation of <Kd>_1, the average spectral attenuation coefficient between surface and first attenuation depth of downwelling planar irradiance
    %In this version of the code, <Kd>_1 is assumed to be known and
    %provided as input in units of [m^-1] In the 2018 paper, it is obtained
    %from a separate neural network algorithm

%% Step 3: Calculation of b, the total scattering coefficient
    %In this version of the code, bp and bw are assumed to be known and
    %provided as input in units of [m^-1] In the 2018 paper, bp is
    %estimated from Chl calculated using OC4v4 model
    b = bp + bw; %[m^-1]

%% Step 4: Calculation of eta, the ratio of the water molecular scattering coefﬁcient to the total scattering coefﬁcient
    eta = bw/b; %[dim]

%% Steps 5 & 7: Calculation of a and bb from Eqs. 9 and 8
    if ~isnan(eta) && ~isnan(muw)
        %find uppermost index of eta and leftmost index of mu values in the
        %LUTs

        %Replacement for LS2_seek_pos that doesn't use a subfunction to
        %find the location of eta in eta look-up table and displays warning
        %if eta value is outside the range of the look-up table. Selected
        %index is the leftmost position of the input look-up table.
        if eta < min(LS2_LUT.eta)
            warning(['eta is outside the lower bound of look-up table.' ...
                ' Solutions of a and bb are output as nan.'])
            a = nan;
            anw = nan;
            bb = nan;
            bbp = nan;
            kappa = nan;
            return
        elseif eta > max(LS2_LUT.eta)
            warning(['eta is outside the upper bound of look-up table.' ...
                ' Solutions of a and bb are output nan.'])
            a = nan;
            anw = nan;
            bb = nan;
            bbp = nan;
            kappa = nan;
            return
        else
            for i = 1:length(LS2_LUT.eta) - 1
                if eta >= LS2_LUT.eta(i) && eta < LS2_LUT.eta(i+1)
                    idx_eta = i;
                end
            end
        end
        
        %Replacement for LS2_seek_pos that doesn't use a subfunction to
        %find the location of muw in muw look-up table and displays warning
        %if muw value is outside the range of the look-up table. Selected
        %indicies is the leftmost position of the input look-up table.
        if muw < min(LS2_LUT.muw)
            warning(['muw is outside the lower bound of look-up table.' ...
                ' Solutions of a and bb are output as nan.'])
            a = nan;
            anw = nan;
            bb = nan;
            bbp = nan;
            kappa = nan;
            return
        elseif muw > max(LS2_LUT.muw)
            warning(['muw is outside the upper bound of look-up table.' ...
                ' Solutions of a and bb are output as nan.'])
            a = nan;
            anw = nan;
            bb = nan;
            bbp = nan;
            kappa = nan;
            return
        else
            for i = 1:length(LS2_LUT.muw) - 1 
                if muw <= LS2_LUT.muw(i) && muw > LS2_LUT.muw(i+1)
                    idx_muw = i;
                end
            end
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
                
        %Calculate using 2-D linear interpolation of a determined
        %from bracketed values of eta and muw
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

        %Calculate bb using 2-D linear interpolation of bb determined
        %from bracketed values of eta and muw
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
    %and recalculate a and bb the same as above Otherwise no correction is
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

            %Calculate using 2-D linear interpolation of a determined
            %from bracketed values of eta and muw
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

            %Calculate bb using 2-D linear interpolation of bb determined
            %from bracketed values of eta and muw
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
    bbw = bw/2;
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
    %Ocean Optics Research Laboratory
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Check function input arguments
    arguments
        bb_a (1,1) double
        lam (1,1) double
        rLUT (101,7) double
    end
    
    %Intepolate minimum and maximum allowable bb/a values that calculate
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