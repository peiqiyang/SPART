%%
% *************************************************************************
% *     functions: Run the SPART model                                    *
% *     Authors:   Peiqi Yang (p.yang@utwente.nl)                         *
% *     Create:    04/Dec/2019                                            *
% *     Update:    12/May/2020                                            *
% *     Faculty of Geo-Information Science and Earth Observation (ITC)    *
% *     University of Twente, 7500 AE Enschede, The Netherlands           *
% *************************************************************************

% Soil-Plant-ATmosphere model (SPART) for top-of-canopy and top-of-atmosphere reflectance 
% Developed by Peiqi Yang (p.yang@utwente.nl),ITC, University of Twente
% Christiaan van der Tol and Wout Verhoef
% Coupling BSM, PROSAIL and SMAC to simulate TOA reflectance

% Updates:
% PY, 11/MAY/2020, optimize I/O. 

%% 0. start fresh and get constants
% dbstop if error
clear 
clc
close all
restoredefaultpath

pathnames                   = 'src_SPART';
addpath(pathnames)          % add the path where the src is and will be removed in the end
%% 1. load the input data and prepare for simulation 

% NOTE: make sure the order of the outputs of the function 'set_input_from_excel' are
% consistent what stated here. 
inputfilename           =   'SPART_input_user_defined.xlsx';
[soilpar,leafbio, canopy, atm, angles, spectral, optipar,sensor]   =   set_input_from_excel(inputfilename);  
% [soilpar,leafbio, canopy, atm,angles, spectral, optipar]  =   set_input_manually();       

nsim                    =   length(leafbio.Cab);
nwl                     =   length(spectral.wlSensor);
[R_TOC,R_TOA,L_TOA]     =   deal(zeros(nwl,nsim));
params_inputs           =   zeros(23,nsim);

%% 2. run the model
tic
for k =1:nsim
[soilpar_i,leafbio_i,canopy_i,atm_i,angles_i]   =      select_parameter(soilpar,leafbio,canopy,atm,angles,k);
[R_TOC(:,k),R_TOA(:,k),L_TOA(:,k)]              =      SPART_main(soilpar_i,leafbio_i,canopy_i,atm_i,angles_i,spectral,optipar);

params_inputs(:,k) = [soilpar_i.B,soilpar_i.lat,soilpar_i.lon,soilpar_i.SMp,...                          % 4 input soil 
    leafbio_i.Cab,leafbio_i.Cdm,leafbio_i.Cw,leafbio_i.Cs,leafbio_i.Cca,leafbio_i.Cant,leafbio_i.N,...  % 7 input leaf
    canopy_i.LAI,canopy_i.LIDFa,canopy_i.LIDFb,canopy_i.hot,...                                         % 4 input canopy
    atm_i.Pa,atm_i.aot550,atm_i.uo3,atm_i.alt_m,atm_i.Pa0...                                            % 5 input atmopshere
    angles_i.tts,angles_i.tto,angles_i.psi];                                                            % 3 input angles

end
toc
output{1} = R_TOC;
output{2} = R_TOA;
output{3} = L_TOA;


%% 3. save output
out_dir_field   =   'testing';
out_dir         =   ['output\',sensor.mission,'-',sensor.name,'-',out_dir_field,'\'];
output_files(out_dir,params_inputs,output,sensor)


%% 4. plot some results 
% the wls are not neccessarily ascending. sorting is needed
[wl,id_sort]   =    sort(spectral.wlSensor);
figure
p1=plot(wl,R_TOC(id_sort,:),'-xk');
hold on
p2= plot(wl,R_TOA(id_sort,:),'-or');
legend([p1(1) p2(1)],'TOC R','TOA R','location','northwest');xlabel('wavelength(nm)'); ylabel('reflectance');
set(gcf, 'units','centimeters','position',[3,3,20,15],'color','w')
title([sensor.mission,'-',sensor.name]);

%% 5. finished, remove the pathes added, to aviod conflicts. 
rmpath(pathnames)
