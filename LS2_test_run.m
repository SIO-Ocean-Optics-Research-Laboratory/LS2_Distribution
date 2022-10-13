%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run LS2 code for ditribution. Utilizes ten inputs and outputs to .xlsx
%file for comparison against xxx.xlsx.

%Created: October 12, 2022
%Completed: N/A
%Updates: N/A
%
%Matthew Kehrli
%Ocean Optics Research Laboratory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clear command window and workspace; close figures
clc; clearvars; close all;

%Define input parameters

%Solar zenith angle input
input_sza = [58.3534804715254; 57.8478623967075; 55.7074826164700; ...
    56.3604406592725; 56.4697386744023; 56.8356539563103; ...
    46.7609493705841; 42.9575051280531; 39.1258070531710; ...
    36.7202617271538];

%Wavlength input
input_lambda = [412, 443, 490, 510, 555, 670];

%Remote sensing reflectance input 10 samples at six wavlengths
input_Rrs = [0.00233524115013000, 0.00294175586230000, ...
    0.00393113794469000, 0.00340969897475000, 0.00229711245356000, ...
    0.000126252474968051;0.00389536285161000,0.00380553119300000, ...
    0.00370193558239000, 0.00287218595770000, 0.00161623731374000, ...
    6.63224220817391e-05; 0.00445894839227000, 0.00429316064858000, ...
    0.00415100665860000, 0.00310168437705000, 0.00160982272779000, ...
    7.69489021946615e-05; 0.00369218043611000, 0.00338770153633000, ...
    0.00335975955251000, 0.00277502852091000, 0.00163977701600000, ...
    9.94781359675408e-05; 0.00439946721488000, 0.00410361567739000, ...
    0.00384560916516000, 0.00281071322991000, 0.00145740677722000, ...
    7.57133335455754e-05; 0.00407536933494000, 0.00393833413124000, ...
    0.00377045558470000, 0.00294465987184000, 0.00155830890220000, ...
    5.98787860868771e-05; 0.00839651952650000, 0.00703053123029000, ...
    0.00515382893062000, 0.00328301126479000, 0.00143446649046000, ...
    7.37164311815052e-05; 0.0108879293669300, 0.00844627966419000, ...
    0.00598548662295000, 0.00380069346088000, 0.00162622388794000, ...
    3.81004114285704e-05; 0.00996484275824000, 0.00805847039179000, ...
    0.00569828863615000, 0.00339472028869000, 0.00146868678397000, ...
    3.54625174513673e-05; 0.00654238127443000, 0.00606473089160000, ...
    0.00511788480497000, 0.00358608311014000, 0.00161083692320000, ...
    3.81378558561512e-05];

%Input of diffuse attenuation coefficients
input_Kd = [0.187632830961948, 0.138869886516291, 0.0903870213746760, ...
    0.0887979802802584, 0.105697824445020, 0.562373761609887; ...
    0.0959924017590809, 0.0776637266248738, 0.0621434486477123, ...
    0.0679142033853613, 0.0929498240015621, 0.541582660606669; ...
    0.0838916051527420, 0.0676778528875712, 0.0548092373676583, ...
    0.0606911584475843, 0.0867712558516902, 0.532867854871956; ...
    0.107015673766332, 0.0859385855948466, 0.0669716298720907, ...
    0.0720650641466051, 0.0953025020320744, 0.535712492953156; ...
    0.0819636283043127, 0.0664145944821315, 0.0545267015934212, ...
    0.0608791413489761, 0.0875697338410355, 0.536433237168806; ...
    0.0870883103721925, 0.0704351866270744, 0.0568638096692833, ...
    0.0627223237995291, 0.0883024973409085, 0.531148295108121; ...
    0.0397738921806664, 0.0355556278806364, 0.0371061601830686, ...
    0.0462176494000053, 0.0758945231309703, 0.500096222235249; ...
    0.0358540683547611, 0.0323961377543049, 0.0353061269500453, ...
    0.0446719301032437, 0.0738646909883413, 0.484404799921858; ...
    0.0352724958202497, 0.0316083929271783, 0.0342524209435186, ...
    0.0435171225162597, 0.0732685592707990, 0.487660368105567; ...
    0.0497417842712074, 0.0416885187358688, 0.0380380139571626, ...
    0.0449932910491427, 0.0717349865479659, 0.474223095982504];

%Input of pure seawater spectral absorption coefficient
input_aw = [0.00467300000000000, 0.00721000000000000, ...
    0.0150000000000000, 0.0325000000000000,0.0592000000000000, ...
    0.439000000000000];

%Input bw
input_bw = [0.00658572000000000, 0.00477700000000000, ...
    0.00309800000000000,0.00259800000000000,0.00184881000000000, ...
    0.000800000000000000];

%Input bp
input_bp = [0.370479484460264, 0.344554283516092, 0.311505199178834, ...
    0.299289309014959, 0.275022608284016, 0.227817235220342; ...
    0.222636201160542, 0.207056692727185, 0.187196152812537, ...
    0.179855127212045, 0.165272279059717, 0.136904649071855; ...
    0.193719314154925, 0.180163335060562, 0.162882362105773, ...
    0.156494818493782, 0.143806049426719, 0.119122921540043; ...
    0.260143039084721, 0.241938898652156, 0.218732514495725, ...
    0.210154768829226, 0.193115192978207, 0.159968555377470; ...
    0.179080613026910, 0.166549012566787, 0.150573903198136, ...
    0.144669044249190, 0.132939121742499, 0.110121212786697; ...
    0.204161647168565, 0.189874940481825, 0.171662446190711, ...
    0.164930585555782, 0.151557835375583, 0.125544177064849; ...
    0.0934854453083219, 0.0869435744176718, 0.0786040887082217, ...
    0.0755215754255463, 0.0693982044450966, 0.0574865723388487; ...
    0.0886492843355160, 0.0824458355445431, 0.0745377656045563, ...
    0.0716147159730051, 0.0658081173805993, 0.0545126942481083; ...
    0.0800160719773182, 0.0744167531707790, 0.0672788197033778, ...
    0.0646404346169708, 0.0593993182966759, 0.0492039129173957; ...
    0.132749681531226, 0.123460200430847, 0.111618099573194, ...
    0.107240919197775, 0.0985457095330905, 0.0816311474490525];

%Input LS2
input_LS2_LUT = load('LS2_LUT.mat');
input_LS2_LUT = input_LS2_LUT.LS2_LUT;

%Input Raman Flag
input_Flag_Raman = 1;

%Preallocate output variables
output_a = nan(size(input_Rrs));
output_anw = nan(size(input_Rrs));
output_bb = nan(size(input_Rrs));
output_bbp = nan(size(input_Rrs));
output_kappa = nan(size(input_Rrs));

%Run LS2 for all samples at every wavelength
for i = 1:size(input_Rrs,1)
    for j = 1:size(input_Rrs,2)
        [output_a(i,j), output_anw(i,j), output_bb(i,j), ...
            output_bbp(i,j), output_kappa(i,j)] = LS2main(input_sza(i), ...
            input_lambda(j), input_Rrs(i,j), input_Kd(i,j), input_aw(j), ...
            input_bw(j), input_bp(i,j), input_LS2_LUT, input_Flag_Raman);
    end
end

%Save inputs and outputs into an excel file
for i = 1:size(input_Rrs,2)
            T = table(repmat(input_lambda(i),size(input_Rrs,1),1), ...
                input_sza, input_Rrs(:,i), input_Kd(:,i), ...
                repmat(input_aw(i),size(input_Rrs,1),1), ...
                repmat(input_bw(i),size(input_Rrs,1),1), ...
                input_bp(:,i),output_a(:,i), output_anw(:,1), ...
                output_bb(:,i), output_bbp(:,i), ...
                output_kappa(:,i));
            T.Properties.VariableNames = {'Input wavelength [nm]' 'Input sza [deg]' 'Input Rrs [1/sr]' ...
                'Input Kd [1/m]' 'Input aw [1/m]' 'Input bw [1/m]' 'Input bp [1/m]' 'Ouput a [1/m]' ...
                'Output anw [1/m]' 'Output bb [1/m]' 'Output bbp [1/m]' 'Output kappa [dim]'};
            outfile=[cd '\LS2_test_run_output'];
            writetable(T,outfile,'FileType','spreadsheet','Sheet',[num2str(input_lambda(i)) ' nm'])
end
