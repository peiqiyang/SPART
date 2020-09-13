function [soilpar,leafbio, canopy, atm,...              % soil, leaf, canopy and atmopsheric parameters
    angles, spectral, optipar, sensor]...         % spectral and angles
    =      set_input_from_excel(inputfilename)

[spectral]          =      define_bands();

%% ===============Section 1 inputs parameters for the models===============
% some parameters are extracted from the excel sheet
[d, N]              = xlsread(inputfilename,'inputparameters'); % read inputparameters
sensorname          = char(N(3,6));                             % the name of the sensor
lut_str             = char(N(4,6));                             % lut yes or no
opt_LUT             = strcmpi(lut_str,'yes');                   % 1 or 0;
str_var             = N(6:end,2);                               % parameter labels

n_val               = sum(~isnan(d),2);                         % number of input values for each parameters (with explaintion lines)
max_nval            = max(n_val);                               % maximal number of input values, size(d,2);
% nsim                = max_nval;                               % number of simulations, (not using LUT)
str_params          = str_var(n_val>0);                         % parameter labels (after removing some lines)
val_params          = d(n_val>0,:);                             % input values of the parameters (uneven lengths)
n_val_params        = n_val(n_val>0);                           % number of input values for each parameters
n_vars              = size(n_val_params,1);                     % number of variables (should be 26)


if opt_LUT          ==  1
    nsim            =   prod(n_val(n_val~=0));                  % number of simulations
    
    [val_cell,val_cell_filled]    =   deal(cell(1,n_vars));     % initialization (input values in cell, even length)
    val_mat_filled                =   zeros(n_vars,nsim);       % initialization (input values in mat)
    
    for jj=1:n_vars
        val_cell{1,jj}     =   val_params(jj,1:n_val_params(jj))'; % convert the data input cell for ndgrid. 
    end
    
    [val_cell_filled{:}]   =   ndgrid(val_cell{:});               % generated a parameter table (key function)
    
    for jj=1:n_vars
        temp_var                =   squeeze(val_cell_filled{jj});
        val_mat_filled(jj,:)    =   temp_var(:)';                %  n_vars x nsim
    end
    
else
    % the number of input values should either be 0 (null rows), 1, or max_nval.
    if  any(ismember(n_val, [0,1,max_nval])==0)
        error('Error. at least two parameters have mulitple inputs with uneven lengths.')
    end
    val_mat_filled = val_params;    
end

% assign input
soilpar     = struct('B',0,'lat',0,'lon',0,'SMp',0);
leafbio     = struct('Cab',0,'Cdm',0,'Cw',0,'Cs',0,'Cca',0,'Cant',0,'N',0);
canopy      = struct('LAI',0,'LIDFa',0,'LIDFb',0,'hot',0);
atm         = struct('Pa',0,'aot550',0,'uo3',0,'uh2o',0,'alt_m',0,'Pa0',0);
angles      = struct('tts',0,'tto',0,'psi',0);
sensor      = struct('FWHM',0,'DOY',0);

struct_input.soilpar    = soilpar;
struct_input.leafbio    = leafbio;
struct_input.canopy     = canopy;
struct_input.atm        = atm;
struct_input.angles     = angles;
struct_input.sensor     = sensor;

struct_name             = fieldnames(struct_input);
nstruct                 = numel(struct_name);

for s_i     =   1:nstruct
    fn          =   fieldnames(struct_input.(struct_name{s_i}));
    nfield      =   numel(fn);
    for f_j=1:nfield
        match_str           =   strcat(struct_name{s_i},'.',fn{f_j});
        match_ind           =   find(contains(str_params,match_str));
        if isempty(match_ind)||isnan(d(1))
            error('cound not find the variable %s in the input file', match_str)
        end
        val             =  val_mat_filled(match_ind(1),:);
        val(isnan(val)) =  val(1);                      % replace nan with d(1)        
        struct_input.(struct_name{s_i}).(fn{f_j})     =   val;
    end
end

soilpar     =    struct_input.soilpar;
leafbio     =    struct_input.leafbio;
canopy      =    struct_input.canopy;
angles      =    struct_input.angles;
atm         =    struct_input.atm;
sensor      =    struct_input.sensor;
DOY         =    sensor.DOY;
sensor.name =    sensorname;

%% ================Section 2 fixed inputs for the models===================
optiparfile              =    'inputdata\PROSPECT_BSM_inputs\Optipar2020_ProspectD_BSM2019.mat';
exterirradfile           =    'inputdata\TOC2TOA_inputs\Extraterrestrial_irradiance.mat';

X                   =   load(['inputdata\TOC2TOA_inputs\sensors_config_SMAC\spart_sensor_info_',sensor.name,'.mat']);
sensor              =   X.sensor;



X                   =   load(optiparfile);
Ea0_                =   load(exterirradfile);
wlSensor            =   sensor.wl_smac;

% Ea0 and Ea_ Extraterrestrial_irradiance
wl_Ea               =   Ea0_.wl_Ea;         % wl
Ea0_                =   Ea0_.Ea;            % W M-2 nm-1
La_                 =   Extraterrestrial_radiance(Ea0_,DOY);

% sensor SRF
atm.La_             =   Spectral_convolution(wl_Ea,La_,sensor); % extraterrestrial radiance W M-2 nm-1 sr-1


%read the 49 coefficients in smac_soefs table
coef                =   sensor.SMAC_coef;  % read coefficients for SMAC for the sensor
atm.coef            =   coef;                  % fitting coefficients of SMAC for the sensor

% the input for PROSPECT and BSM
optipar             =   X.optipar;              % leaf pigments absorption coefficients


% other inputs, not important for reflectance, but necessary to run the model
leafbio.rho_thermal         = 0.01;             % reflectance in the thermal range
leafbio.tau_thermal         = 0.01;             % transmittance in the thermal range
spectral.wlSensor           = wlSensor';



%%
% struct_input.soilpar    = soilpar;
% struct_input.leafbio    = leafbio;
% struct_input.canopy     = canopy;
% struct_input.angles     = angles;
% struct_input.atm        = atm;
% struct_input.sensor     = sensor;
% struct_input.optipar    = optipar;
% struct_input.spectral   = spectral;

%% using this code to check the spectral convolution results
% figure
% plot(wl_Ea,La_)
% hold on
% plot(wlSensor,atm.La_,'x-')
% ylabel('Extraterrestrial radiance (W M-2 nm-1 sr-1)')
% xlabel('wavelength (nm)')
% legend('before convolution','after convolution')
% xlim([300,1200])

end
