clear
close all
clc
% Soil-Plant-Atmosphere radiative transfer model (SPART) for top-of-canopy and top-of-atmosphere reflectance
% Developed by Peiqi Yang (p.yang@utwente.nl),ITC, University of Twente

% edited by Peiqi Yang
% Coupling BSM, PROSAIL and SMAC to simulate TOA reflectance

% This is an example for local sensitivity analysis of SPART to LAI


%% 1. load the input data
% NOTE: make sure the order of the outputs of the function 'input_data' are
% consistent what stated here. 

inputfilename           =   'SPART_input_user_defined.xlsx';
[soilpar,leafbio, canopy, atm, angles, spectral, optipar,sensor]   =   set_input_from_excel(inputfilename);  

%% Atmosphere inputs
LAI                     =   [0.5,1,3,6];
nsim                    =   length(LAI);
nwl                     =   length(spectral.wlSensor);
[R_TOC,R_TOA,L_TOA]     =   deal(zeros(nwl,nsim));
params_inputs            =   zeros(23,nsim);

tic
for k =1:nsim
[soilpar_i,leafbio_i,canopy_i,atm_i,angles_i]   =      select_parameter(soilpar,leafbio,canopy,atm,angles,1);
canopy_i.LAI                                    =      LAI(k);
[R_TOC(:,k),R_TOA(:,k),L_TOA(:,k)]              =      SPART_main(soilpar_i,leafbio_i,canopy_i,atm_i,angles_i,spectral,optipar);

params_inputs(:,k) = [soilpar_i.B,soilpar_i.lat,soilpar_i.lon,soilpar_i.SMp,...                          % 4 input soil 
    leafbio_i.Cab,leafbio_i.Cdm,leafbio_i.Cw,leafbio_i.Cs,leafbio_i.Cca,leafbio_i.Cant,leafbio_i.N,...  % 7 input leaf
    canopy_i.LAI,canopy_i.LIDFa,canopy_i.LIDFb,canopy_i.hot,...                                         % 4 input canopy
    atm_i.Pa,atm_i.aot550,atm_i.uo3,atm_i.alt_m,atm_i.Pa0...                                            % 5 input atmopshere
    angles_i.tts,angles_i.tto,angles_i.psi];   

end
toc

output{1} = R_TOC;
output{2} = R_TOA;
output{3} = L_TOA;


%% 3. save output
out_dir_field   =   'LSA_LAI';
out_dir         =   ['output\',sensor.mission,'-',sensor.name,'-',out_dir_field,'\'];
output_files(out_dir,params_inputs,output,sensor)

%% 4. make plots
plot_LSA_LAI(out_dir)
