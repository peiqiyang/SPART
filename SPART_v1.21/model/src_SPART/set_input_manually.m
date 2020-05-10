function [soilpar,leafbio, canopy, atm,...      % soil, leaf, canopy and atmopsheric parameters
          angles, spectral, optipar]...         % spectral and angles 
                     =      set_input_manually()   
                 

 [spectral]          =      define_bands(); 
      
%% ===============Section 1 inputs parameters for the models===============
% simulation.cols  = 1;                               % simulation method. 0 or 1
sensorname                          =    'S3A_OLCI';
default_optipar_file                =    'inputdata/SLC_inputs/Optipar2017_ProspectD_BSM.mat';
default_exterirrad_file             =    'inputdata/SLC_inputs/Extraterrestrial_irradiance.mat';
default_wl_file                     =    'inputdata/Sensors/wl_S3.dat';
default_SRF_file                    =    'inputdata/Sensors/sensors_SRF.xlsx';
default_SMAC_coef_file              =    'inputdata/SMAC_COEFS/SMAC_Sentinel3/coef_SENTINEL3A_OLCI_V32_CONTINENTAL.dat';

soilpar.B       = 0.5;
soilpar.lat     = 25;
soilpar.lon     = 45;
soilpar.SMp     = 15;

leafbio.Cab     = 40;               % chlorophyll content               [ug cm-2]
leafbio.Cdm     = 0.01;             % dry matter content                [g cm-2]
leafbio.Cw      = 0.02;             % leaf water thickness equivalent   [cm]
leafbio.Cs      = 0;                % senescent material                [fraction]
leafbio.Cca     = 10;               % carotenoids                       [ug cm-2]
leafbio.Cant    = 10;               % anthocyanins                      [?]
leafbio.N       = 1.5;              % leaf structure parameter (affects the ratio of refl: transmittance) []

canopy.LAI      = 3;        
canopy.LIDFa    = -0.35;
canopy.LIDFb    = -0.15;
canopy.hot      = 0.05;             % ratio between leaf width: canopy height. This is the hot spot parameter. 
                                    % you were probably avoiding the hotspot,
                                    % so in that case the model is insensitive
                                    % to canopy.hot   
                                  
angles.tts      = 40;               % solar zenith angle (deg)
angles.tto      = 0;                % observation zenith angle (deg)
angles.psi      = 0;                % absolute value (deg) of the difference between solar and zenith azimuth angle (arbitrary if angles.tto = 0)

atm.Pa          = -999;             % Air pressure (hPa)
atm.aot550      = 0.3;              % AOT at 550 nm
atm.uo3         = 0.3;              % Ozone content (cm)  0.3 cm= 300 Dobson Units
atm.uh2o        = 3;                % Water vapour (g/cm2)
atm.alt_m       = 0;                % altitude (m). In case, no air pressure input avabile
atm.Pa0         = 1013.25;          % sea level air pressure. In case, no air pressure input avabile

sensor.FWHM     = 3;                % sensor FWHM
DOY             = 100;              % day of year

%% ================Section 2 fixed inputs for the models===================

optipar_filename    =   default_optipar_file;
Ea0_filename        =   default_exterirrad_file;
wl_filename         =   default_wl_file;
nom_smac            =   default_SMAC_coef_file;
SRF_filename        =   default_SRF_file;


X                   =   load(optipar_filename);
wlSensor            =   load(wl_filename);
Ea0_                =   load(Ea0_filename);

% Ea0 and Ea_ Extraterrestrial_irradiance
wl_Ea               =   Ea0_.wl_Ea;         % wl
Ea0_                =   Ea0_.Ea;            % W M-2 nm-1
La_                 =   Extraterrestrial_radiance(Ea0_,DOY,angles.tts);

% sensor SRF
sensor.name         =   sensorname;                             % the name of the sensor
sensor.wl           =   wlSensor';                              % center wavelengths for the sensor
sensor.save         =   1;                                      % do you want to save the SRF created (if it is not a new sensor)
sensor              =   read_or_create_srf(sensor,SRF_filename);   % read or create the spectral response function curve
atm.La_             =   Spectral_convolution(wl_Ea,La_,sensor); % extraterrestrial radiance W M-2 nm-1 sr-1


%read the 49 coefficients in smac_soefs table         
coef                =   read_coeffs(nom_smac);  % read coefficients for SMAC for the sensor
atm.coef            =   coef;                   % fitting coefficients of SMAC for the sensor

% the input for PROSPECT and BSM
optipar            = X.optipar;        % leaf pigments absorption coefficients

% other inputs, not important for reflectance, but necessary to run the model
leafbio.rho_thermal         = 0.01;             % reflectance in the thermal range
leafbio.tau_thermal         = 0.01;             % transmittance in the thermal range
spectral.wlSensor           = wlSensor'; 



end

%% using this code to check the spectral convolution
% figure
% plot(wl_Ea,La_)
% hold on
% plot(wlSensor,atm.La_,'x-')
% ylabel('Extraterrestrial radiance (W M-2 nm-1 sr-1)')
% xlabel('wavelength (nm)')
% legend('before convolution','after convolution')
% xlim([300,1200])