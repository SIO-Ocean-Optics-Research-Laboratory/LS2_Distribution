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
%       to construct the a and bb LUTs LUT.eta [21x1 double]: 21 values of
%       eta (in ascending order) used to construct the a and bb LUTs
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
    %In this version of the code, <Kd>_1 is assumed to be known and provided as input in units of [m^-1]
    %In the 2018 paper, it is obtained from a separate neural network algorithm

%% Step 3: Calculation of b, the total scattering coefficient
    %In this version of the code, bp and bw are assumed to be known and provided as input in units of [m^-1]
    %In the 2018 paper, bp is estimated from Chl calculated using OC4v4 model
    b = bp + bw; %[m^-1]

%% Step 4: Calculation of eta, the ratio of the water molecular scattering coefﬁcient to the total scattering coefﬁcient
    eta = bw/b; %[dim]

%% Steps 5 & 7: Calculation of a and bb from Eqs. 9 and 8
    if ~isnan(eta) && ~isnan(muw)
        
        %find uppermost index of eta and leftmost index of mu values in the LUTs
        %Calls subfunction LS2_seek_pos
        %RR - COULD DO THIS MATLAB - SHOULD ADD A CHECK TO MAKE SURE VALUES OF ETA AND MUW ARE IN RANGE OF TABLE
        idx_eta = LS2_seek_pos(eta,LS2_LUT.eta,'eta');
        idx_muw = LS2_seek_pos(muw,LS2_LUT.muw,'muw');
        
        %calculation of a from Eq. 9
            %a is first calculated for four combinations of LUT values bracketing the lower and upper limits of eta and muw
              a00 = Kd./(LS2_LUT.a(idx_eta,idx_muw,1) + LS2_LUT.a(idx_eta,idx_muw,2).*Rrs + ...
                    LS2_LUT.a(idx_eta,idx_muw,3).*Rrs.^2 + LS2_LUT.a(idx_eta,idx_muw,4).*Rrs.^3);

              a01 = Kd./(LS2_LUT.a(idx_eta,idx_muw+1,1) + LS2_LUT.a(idx_eta,idx_muw+1,2).*Rrs + ...
                    LS2_LUT.a(idx_eta,idx_muw+1,3).*Rrs.^2 + LS2_LUT.a(idx_eta,idx_muw+1,4).*Rrs.^3);

              a10 = Kd./(LS2_LUT.a(idx_eta+1,idx_muw,1) + LS2_LUT.a(idx_eta+1,idx_muw,2).*Rrs + ...
                    LS2_LUT.a(idx_eta+1,idx_muw,3).*Rrs.^2 + LS2_LUT.a(idx_eta+1,idx_muw,4).*Rrs.^3);

              a11 = Kd./(LS2_LUT.a(idx_eta+1,idx_muw+1,1) + LS2_LUT.a(idx_eta+1,idx_muw+1,2).*Rrs + ...
                    LS2_LUT.a(idx_eta+1,idx_muw+1,3).*Rrs.^2 + LS2_LUT.a(idx_eta+1,idx_muw+1,4).*Rrs.^3);
            
            %calculated a is then linearly interpolated between the two values of eta for each value of muw
              tmp1 = LS2_interp_line(LS2_LUT.eta(idx_eta),LS2_LUT.eta(idx_eta+1),a00,a10,eta);
              tmp2 = LS2_interp_line(LS2_LUT.eta(idx_eta),LS2_LUT.eta(idx_eta+1),a01,a11,eta);

              %RR - I THINK THIS BUILT-IN MATLAB FUNCTION COULD BE SUBSTITUTED FOR PERFORMING LINEAR INTERPOLATION
              %tmp1 = interp1(LS2_LUT.eta(idx_eta:idx_eta+1),[a00 a10],eta);
              %tmp2 = interp1(LS2_LUT.eta(idx_eta:idx_eta+1),[a01 a11],eta);
            
            %final a is calculated from linear interpolation between the two values of a for each muw
              a = LS2_interp_line(LS2_LUT.muw(idx_muw),LS2_LUT.muw(idx_muw+1),tmp1,tmp2,muw); %[m^-1]
              %RR a = interp1((LS2_LUT.muw(idx_muw:idx_muw+1),[tmp1 tmp2],muw)); %[m^-1]

        %calculation of bb from Eq. 8
            %bb is first calculated for four combinations of LUT values bracketing the lower and upper limits of eta and muw
            bb00 = Kd.*(LS2_LUT.bb(idx_eta,idx_muw,1).*Rrs + LS2_LUT.bb(idx_eta,idx_muw,2).*Rrs.^2 + ...
                   LS2_LUT.bb(idx_eta,idx_muw,3).*Rrs.^3);

            bb01 = Kd.*(LS2_LUT.bb(idx_eta,idx_muw+1,1).*Rrs + LS2_LUT.bb(idx_eta,idx_muw+1,2).*Rrs.^2 + ...
                   LS2_LUT.bb(idx_eta,idx_muw+1,3).*Rrs.^3);

            bb10 = Kd.*(LS2_LUT.bb(idx_eta+1,idx_muw,1).*Rrs + LS2_LUT.bb(idx_eta+1,idx_muw,2).*Rrs.^2 + ...
                   LS2_LUT.bb(idx_eta+1,idx_muw,3).*Rrs.^3);

            bb11 = Kd.*(LS2_LUT.bb(idx_eta+1,idx_muw+1,1).*Rrs + LS2_LUT.bb(idx_eta+1,idx_muw+1,2).*Rrs.^2 + ...
                   LS2_LUT.bb(idx_eta+1,idx_muw+1,3).*Rrs.^3);

        %bb is then linearly interpolated between the two values of eta for different muw
            tmp1 = LS2_interp_line(LS2_LUT.eta(idx_eta),LS2_LUT.eta(idx_eta+1),bb00,bb10,eta);
            tmp2 = LS2_interp_line(LS2_LUT.eta(idx_eta),LS2_LUT.eta(idx_eta+1),bb01,bb11,eta);
            
            %RR
            %tmp1 = interp1(LS2_LUT.eta(idx_eta:idx_eta+1),[bb00 bb10],eta);
            %tmp2 = interp1(LS2_LUT.eta(idx_eta:idx_eta+1),[bb01 bb11],eta);

        %final bb is calculated from linear interpolation between the two values of muw    
            bb = LS2_interp_line(LS2_LUT.muw(idx_muw),LS2_LUT.muw(idx_muw+1),tmp1,tmp2,muw); %[m^-1]

            %RR bb = interp1((LS2_LUT.muw(idx_muw:idx_muw+1),[tmp1 tmp2],muw)); %[m^-1]

     else
        a = NaN; bb = NaN; %return NaNs if either eta or muw are NaNs
    end

%% Step 9: Application of the Raman scattering correction if selected
    %If Flag_Raman is set to 1 (true), apply Raman correction to input Rrs and recalculate a and bb the same as above
    %Otherwise no correction is applied and original values are returned with kappa value of 1
    if Flag_Raman 
        kappa = LS2_calc_kappa(bb/a,lambda,LS2_LUT.kappa); %Calls subfunction LS2_calc_kappa

        %Apply Raman scattering correction to Rrs and recalculate a & bb
        if ~isnan(kappa)
            Rrs = Rrs.*kappa;
            
            %calculation of a from Eq. 9
                %a is first calculated for four combinations of LUT values bracketing the lower and upper limits of eta and muw
                a00 = Kd./(LS2_LUT.a(idx_eta,idx_muw,1) + LS2_LUT.a(idx_eta,idx_muw,2).*Rrs + ...
                    LS2_LUT.a(idx_eta,idx_muw,3).*Rrs.^2 + LS2_LUT.a(idx_eta,idx_muw,4).*Rrs.^3);

                a01 = Kd./(LS2_LUT.a(idx_eta,idx_muw+1,1) + LS2_LUT.a(idx_eta,idx_muw+1,2).*Rrs + ...
                    LS2_LUT.a(idx_eta,idx_muw+1,3).*Rrs.^2 + LS2_LUT.a(idx_eta,idx_muw+1,4).*Rrs.^3);

                a10 = Kd./(LS2_LUT.a(idx_eta+1,idx_muw,1) + LS2_LUT.a(idx_eta+1,idx_muw,2).*Rrs + ...
                    LS2_LUT.a(idx_eta+1,idx_muw,3).*Rrs.^2 + LS2_LUT.a(idx_eta+1,idx_muw,4).*Rrs.^3);

                a11 = Kd./(LS2_LUT.a(idx_eta+1,idx_muw+1,1) + LS2_LUT.a(idx_eta+1,idx_muw+1,2).*Rrs + ...
                    LS2_LUT.a(idx_eta+1,idx_muw+1,3).*Rrs.^2 + LS2_LUT.a(idx_eta+1,idx_muw+1,4).*Rrs.^3);
            
            %calculated a is then linearly interpolated between the two values of eta for each value of muw
                tmp1 = LS2_interp_line(LS2_LUT.eta(idx_eta),LS2_LUT.eta(idx_eta+1),a00,a10,eta);
                tmp2 = LS2_interp_line(LS2_LUT.eta(idx_eta),LS2_LUT.eta(idx_eta+1),a01,a11,eta);

                %RR - I THINK BUILT-IN MATLAB FUNCTION COULD BE SUBSTITUTED FOR PERFORMING LINEAR INTERPOLATION
                %tmp1 = interp1(LS2_LUT.eta(idx_eta:idx_eta+1),[a00 a10],eta);
                %tmp2 = interp1(LS2_LUT.eta(idx_eta:idx_eta+1),[a01 a11],eta);
            
            %final a is calculated from linear interpolation between the two values of a for each muw
                 a = LS2_interp_line(LS2_LUT.muw(idx_muw),LS2_LUT.muw(idx_muw+1),tmp1,tmp2,muw); %[m^-1]
                %RR a = interp1((LS2_LUT.muw(idx_muw:idx_muw+1),[tmp1 tmp2],muw)); %[m^-1]

            %calculation of bb from Eq. 8
            %bb is first calculated for four combinations of LUT values bracketing the lower and upper limits of eta and muw
                bb00 = Kd.*(LS2_LUT.bb(idx_eta,idx_muw,1).*Rrs + LS2_LUT.bb(idx_eta,idx_muw,2).*Rrs.^2 + ...
                   LS2_LUT.bb(idx_eta,idx_muw,3).*Rrs.^3);

                bb01 = Kd.*(LS2_LUT.bb(idx_eta,idx_muw+1,1).*Rrs + LS2_LUT.bb(idx_eta,idx_muw+1,2).*Rrs.^2 + ...
                   LS2_LUT.bb(idx_eta,idx_muw+1,3).*Rrs.^3);

                bb10 = Kd.*(LS2_LUT.bb(idx_eta+1,idx_muw,1).*Rrs + LS2_LUT.bb(idx_eta+1,idx_muw,2).*Rrs.^2 + ...
                   LS2_LUT.bb(idx_eta+1,idx_muw,3).*Rrs.^3);

                bb11 = Kd.*(LS2_LUT.bb(idx_eta+1,idx_muw+1,1).*Rrs + LS2_LUT.bb(idx_eta+1,idx_muw+1,2).*Rrs.^2 + ...
                   LS2_LUT.bb(idx_eta+1,idx_muw+1,3).*Rrs.^3);

            %bb is then linearly interpolated between the two values of eta for different muw
                tmp1 = LS2_interp_line(LS2_LUT.eta(idx_eta),LS2_LUT.eta(idx_eta+1),bb00,bb10,eta);
                tmp2 = LS2_interp_line(LS2_LUT.eta(idx_eta),LS2_LUT.eta(idx_eta+1),bb01,bb11,eta);
            
                %RR
                %tmp1 = interp1(LS2_LUT.eta(idx_eta:idx_eta+1),[bb00 bb10],eta);
                %tmp2 = interp1(LS2_LUT.eta(idx_eta:idx_eta+1),[bb01 bb11],eta);

            %final bb is calculated from linear interpolation between the two values of muw    
                bb = LS2_interp_line(LS2_LUT.muw(idx_muw),LS2_LUT.muw(idx_muw+1),tmp1,tmp2,muw); %[m^-1]
                %RR bb = interp1((LS2_LUT.muw(idx_muw:idx_muw+1),[tmp1 tmp2],muw)); %[m^-1]
        end

    else %if Flag is not 1, do nothing and return original values of a, anw, bb, bbp with kappa returned as 1
        kappa = 1;
    
    end

%% Steps 6 & 8: Calculation of anw and bbp
    anw = a - aw; % spectral nonwater absorption coefficient [m^-1]
    bbw = bw/2;
    bbp = bb - bbw; %spectral particulate backscattering coefficient [m^-1]

%% If output coefficients are negative, replace with NaN
    a(a<0) = NaN;   anw(anw<0)=NaN;
    bb(bb<0) = NaN; bbp(bbp<0)=NaN;

end
%end of main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Additional subfunctions that are called
function idx = LS2_seek_pos(param,lut,type)
    %MK remake of LS2 seek_pos subroutine. The subroutine finds the leftmost
    %position of the input parameter in relation to its input LUT.
    %
    %Inputs: param, lut, type
    %   param (1x1 Double): Input muw or eta value
    %
    %   lut (nx1 Double): Look-up table of muw or eta values used to determine
    %   coefficients in Loisel et al. 2018. If the input is associated with muw
    %   the lut must be 8x1 and sorted in descending order, and if the input is
    %   associated with eta the lut must be 21x1 and sorted in ascending order.
    %
    %   type (String): Characterize param input. Valid values are 'muw' or
    %   'eta'. Other inputs will produce and error.
    %
    %Output: idx
    %   idx (1x1 Double): Leftmost position/index of input param in relation to
    %   its LUT.
    %
    %Created: September 8, 2021
    %Completed: September 8, 2021
    %Updates: 
    %
    %Matthew Kehrli
    %Ocean Optics Research Laboratory
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Check function input arguments
    arguments 
        param (1,1) double
        lut (:,1) double
        type (1,1) string
    end
    
    %mu_w lookup
    if strcmp(type,'muw')
        if length(lut) ~= 8 || ~issorted(lut,'descend')
            error(['Look-up table for mu_w must be a 8x1 array sorted in '...
                'descending order'])
        end
          
       if param >= lut(1)
           idx = 1;
       end
       if param <=lut(end)
           idx = length(lut) - 1;
       end
       for i = 1:length(lut) - 1
           if param <= lut(i) && param > lut(i+1)
               idx = i;
           end
       end
       
    %Eta lookup
    elseif strcmp(type,'eta')
        if length(lut) ~= 21 || ~issorted(lut)
            error(['Look-up table for eta must be a 21x1 array sorted in '...
                'ascending order'])
        end
        if param <= lut(1)
            idx = 1;
        end
        if param >= lut(end)
            idx = length(lut) - 1;
        end
        for i = 1:length(lut) - 1
            if param >= lut(i) && param < lut(i+1)
                idx = i;
            end
        end
        
    end
end

function yq = LS2_interp_line(x1,x2,y1,y2,xq)
    %MK remake of LS2 interp_line subroutine. The subroutine linearly
    %interpolates/extrapolates a value of yq for a given query point, xq.
    %
    %Inputs: x1,x2,y1,y2,xq
    %   x1 (1x1 Double): First x sample point.
    %
    %   x2 (1x1 Double): Second x sample point.
    %
    %   y1 (1x1 Double): First y sample point.
    %
    %   y2 (1x1 Double): Second y sample point.
    %
    %   xq (1x1 Double): Query x point.
    %
    %Output: yq
    %   yq (1x1 Double): Query y point.
    %
    %Created: September 8, 2021
    %Completed: September 8, 2021
    %Updates: 
    %
    %Matthew Kehrli
    %Ocean Optics Research Laboratory
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Check function input arguments
    arguments 
        x1 (1,1) double
        x2 (1,1) double
        y1 (1,1) double
        y2(1,1) double
        xq (1,1) double
    end
    
    % Find piecewise slopes
    a = (y2 - y1)./(x2 - x1);
    
    % Find piecewise intercepts
    b = y1 - (a.*x1);
    
    % Find functional values at query points
    yq = a*xq + b;
end

function kappa = LS2_calc_kappa(bb_a,lam,rLUT)
    %MK remake of LS2 calc_kappa subroutine. The subroutine determines kappa
    %using a linear interpolation/extrapolation from the Raman scattering
    %look-up tables.
    %
    %Inputs: Rrs, Kd, eta, muw, pos_eta, pos_mu, lut_eta, lut_mu, lut_a
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Check function input arguments
    arguments
        bb_a (1,1) double
        lam (1,1) double
        rLUT (101,7) double
    end
    
    %Intepolate minimum and maximum allowable bb/a values that calculate kappa
    %to the input wavelegnth.
    mins = interp1(rLUT(:,1),rLUT(:,6),lam,'linear');
    maxs = interp1(rLUT(:,1),rLUT(:,7),lam,'linear');
    
    %Check if bb/a ratio falls within min/max range
    if bb_a >= mins && bb_a <= maxs
        %Calculate all kappas for the input bb/a ratio
        kappas = rLUT(:,2).*bb_a.^3 + rLUT(:,3).*bb_a.^2 + rLUT(:,4).*bb_a ...
            + rLUT(:,5);
        %Calculate output kappa using a 1-D linear interpolation to the input
        %wavlength
        kappa = interp1(rLUT(:,1),kappas,lam,'linear');
    else
        kappa = nan;
        %Send warning message to user that no Raman correction is applied to
        %input
        warning(['No Raman Correction since bb/a value is outside of the' ...
            ' acceptable range. Kappa set to nan and no correction is' ...
            ' applied. See Raman Correction LUT.'])
    end
end
