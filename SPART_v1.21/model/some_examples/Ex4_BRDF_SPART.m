clear
close all
clc
% Soil-Plant-Atmosphere radiative transfer model (SPART) for top-of-canopy and top-of-atmosphere reflectance
% Developed by Peiqi Yang (p.yang@utwente.nl),ITC, University of Twente

% edited by Peiqi Yang
% Coupling BSM, PROSAIL and SMAC to simulate TOA reflectance

% This is an example for BRDF simulation



%% 1. load the input data
% NOTE: make sure the order of the outputs of the function 'input_data' are
% consistent what stated here. 


inputfilename           =   'SPART_input_user_defined.xlsx';
[soilpar,leafbio, canopy, atm, angles, spectral, optipar,sensor]   =   set_input_from_excel(inputfilename); 

%% BRDF inputs
tts                 = angles.tts(1);
anglesfile          = load('inputdata/BRDF_files/brdf_angles2.dat'); %Multiple observation angles in case of BRDF calculation
directional.tto     = anglesfile(:,1);              % [deg]              %Observation zenith Angles for calcbrdf
directional.psi     = anglesfile(:,2);              % [deg]              %Observation zenith Angles for calcbrdf
directional.noa     = length(directional.tto);      %                    %Number of Observation Angles


psi_hoversampling   = [0    ; 0     ;0     ;0     ;0    ;2  ;358];  % [noa_o]           angles for hotspot oversampling
tto_hoversampling   = [tts  ; tts+02;tts+04;tts-02;tts-4;tts;tts];  % [noa_o]           angles for hotspot oversampling
noah_o              = size(tto_hoversampling,1);                    % [1]               number of oversampling angles
psi_poversampling   = [000*ones(6,1);180*ones(6,1);090*ones(6,1);270*ones(6,1)];%       angles for plane oversampling
tto_poversampling   = [10:10:60     , 10:10:60    , 10:10:60    , 10:10:60]';   %       angles for plane oversampling
noap_o              = size(tto_poversampling,1);                    % [1]               number of oversampling angles


directional.psi     = [directional.psi;psi_hoversampling;psi_poversampling];   % [..]   observer azimuth angle
directional.tto     = [directional.tto;tto_hoversampling;tto_poversampling];   % [..]   observer zenith  angle
directional.noa     = directional.noa+noap_o+noah_o;


tic
nsim                    =   directional.noa;
nwl                     =   length(spectral.wlSensor);
[R_TOC,R_TOA,L_TOA]     =   deal(zeros(nwl,nsim));
parms_inputs            =   zeros(23,nsim);

for k =1:nsim 
    
[soilpar_i,leafbio_i,canopy_i,atm_i,angles_i]    =      select_parameter(soilpar,leafbio,canopy,atm,angles,1);
    angles_i.tto                                 =      directional.tto(k);
    angles_i.psi                                 =      directional.psi(k);
[R_TOC(:,k),R_TOA(:,k),L_TOA(:,k)]               =      SPART_main(soilpar_i,leafbio_i,canopy_i,atm_i,angles_i,spectral,optipar);
    
    
parms_inputs(:,k) = [soilpar_i.B,soilpar_i.lat,soilpar_i.lon,soilpar_i.SMp,...                          % 4 input soil 
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
out_dir_field   =   'BRDF';
out_dir         =   ['output\',sensor.mission,'-',sensor.name,'-',out_dir_field,'\'];
output_files(out_dir,parms_inputs,output,sensor)


tto = directional.tto;
psi = directional.psi;
wl  = sensor.wl_smac;
band_2_plot = 3;
plot_BRDF_function(tto,psi,R_TOC,R_TOA,wl,3)

