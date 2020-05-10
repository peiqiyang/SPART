function [R_TOC,R_TOA,L_TOA] = SPART_main(soilpar,leafbio,canopy,atm,angles,spectral,optipar)
%% Soil-Plant-Atmosphere Radiative Transfer model for top-of-canopy and top-of-atmosphere reflectance 
% Developed by Peiqi Yang               (p.yang@utwente.nl)
%              Christiaan van der Tol   (c.vandertol@utwente.nl)
%              Wout Verhoef             (w.verhoef@utwente.nl)


% University of Twente, Faculty of Geo-Information Science and Earth Observation (ITC), 
% Department of Water Resources, 
% P.O. Box 217, 7500 AE Enschede, The Netherlands
% Corresponding author: peiqi Yang
% +31(0)685771997;
% p.yang@utwente.nl,peiqiyangweb@gmail.com

%% main inputs
% soil properties: 
% soilpar.B         =   soil brightness [0,1] unitless
% soilpar.lat       =   soil spectral shape parameters [25-45] unitless
% soilpar.lon       =   soil spectral shape parameters [30-60] unitless
% soilpar.SMp       =   soil mositure in percentage [5-55] unitless

% leaf (bio)physical properties:
% leafbio.Cab       =   Chlorophyll content [0- 80], ug cm-2 
% leafbio.Cdm       =   leaf mass per unit area [0-0.02] g cm-2
% leafbio.Cw        =   Equivalent water thickness [0-0.1] cm
% leafbio.Cca       =   Caratenoid content [0-30] ug cm-2
% leafbio.Cs        =   brown pigments [0-1] unitless, 0=green 1=brown
% leafbio.Cant      =   Anthocyanin content [0-30] ug cm-2
% leafbio.N         =   leaf structure parameter [0-3] related to thickness

% canopy structural properties
% canopy.LAI        =   Leaf area index [0,7]  m m-2
% canopy.LIDFa      =   leaf inclination parameter a [-1,1]
% canopy.LIDFb      =   parameter b [-1,1], abs(a)+abs(b)<=1
% canopy.hot        =   hot spot parameters, [0-0.2] leaf width/canopy height

% atmopsheric properties
% atm.Pa            =   Air pressure [500-1300] hPa/estimated from elevation
% atm.aot550        =   aerosol optical thickness at 550nm [0-2] unitless
% atm.uo3           =   Ozone content [0-0.8] cm-atm or [0-0.0171] kg m-2
% atm.uh2o          =   Water vapour  [0-8.5] g cm-2 or [0-85] kg m-2
% atm_alt_m         =   groudn elevation in meters in case no Pa available
% atm.Ea            =   extraterrestrial irradiance spectral

% sun-observer geometry
% angles.tts        =   solar zenith angle
% angles.tto        =   viewing zenith angle
% angles.psi        =   relative azumith angle between the sun and viewer

%% main ouputs
% R_TOC     Top-of-Canopy reflectance
% R_TOA     Top-of-Atmopshere reflectance
% L_TOA     Top-of-Atmopshere radiance


% main functions 
% - BSM (Brightness-shape-mositure)
% - PROSAIL
% - SAIL (SAILH, with hotspot effects)
% - SMAC (modified SMAC model)

% Coupling BSM, PROSAIL and SMAC to simulate TOA reflectance 


% For every complex question there is a simple and wrong solution

% Updates: -created                                     29/Aug/2019 by P.Y.
% Updates: -add spectral convolution                    30/Aug/2019 by P.Y.
% Updates: -I/O improvement                             2/Sept/2019 by P.Y.
% Updates: -SAIL and SMAC interaction                   22/Sept/2019 by P.Y. 
% Updates: -SAIL and SMAC interaction improvement       07/Oct/2019 by P.Y.
%           1. tg was not seperated into tgs and tgv (two-ways gasous transmitamission into one by one)
%           2. SMAC modification, to get tss and too excluding the effects of
%               gas, such as O3 and H2O
%           3. rso was re-caculated                     17/Oct/2019 by P.Y.
%Updates:  -compute TOA radiance                        16/Nov/2019 by P.Y.

%%
%% 0. get prepared
% spectral information and optical parameters for PROSPECT and BSM
% optipar         =  other_input{1};
% spectral        =  other_input{2};

% the fixed optical parameters for BSM
soilemp.SMC     =  25;               % soil moisture content
soilemp.film    =  0.015;            % water film optical thickness

soilspec.GSV    =  optipar.GSV;      % Global vectors for dry soil reflectance
soilspec.kw     =  optipar.Kw;       % water absorption coefficients
soilspec.nw     =  optipar.nw;       % refractive index                

% spectral information 
IwlP            =   spectral.IwlP;
IwlT            =   spectral.IwlT;
nwlP            =   spectral.nwlP;
nwlT            =   spectral.nwlT;

%% 1 run the soil reflectance model
[rho,tau,rs]    =   deal(zeros(nwlP + nwlT,1));
rs(IwlP)        =   BSM(soilpar,soilspec,soilemp);
rs(IwlT)        =   rs(nwlP) * ones(nwlT,1);
tau(IwlT)       =   leafbio.tau_thermal;
rho(IwlT)       =   leafbio.rho_thermal;
soil.refl       =   rs;

%%
% soilR(1:2000)= rs(1:2000);
% soilR(2001:2101) = rs(2001);
% wl = 400:2500;
%%
% keyboard
%% 2 run the leaf radiative transfer model
[leafopt]       =   PROSPECT_5D(leafbio,optipar);
rho(IwlP)       =   leafopt.refl;
tau(IwlP)       =   leafopt.tran;
leafopt.refl    =   rho;                  
leafopt.tran    =   tau;

%% 3 run the canopy radiative transfer model (TOC reflectance factors)
rad             =   SAILH(soil,leafopt,canopy,angles);
rv_so           =   interp1(spectral.wlS, rad.rso ,spectral.wlSensor,'splines',1E-4);
rv_do           =   interp1(spectral.wlS, rad.rdo ,spectral.wlSensor,'splines',1E-4);
rv_dd           =   interp1(spectral.wlS, rad.rdd ,spectral.wlSensor,'splines',1E-4);
rv_sd           =   interp1(spectral.wlS, rad.rsd ,spectral.wlSensor,'splines',1E-4);


%% 4 run the atmopshere radiative transfer model
atmopt          = SMAC(angles,atm);
ta_ss           = atmopt.Ta_ss;
ta_sd           = atmopt.Ta_sd; 
ta_oo           = atmopt.Ta_oo;
ta_do           = atmopt.Ta_do;
ra_dd           = atmopt.Ra_dd;
ra_so           = atmopt.Ra_so;
T_g             = atmopt.Tg;

%% 5. upscale TOC to TOA
R_TOC           =   (ta_ss.*rv_so+ta_sd.*rv_do)./(ta_ss+ta_sd);
rtoa1           =   (ta_sd.*rv_do+ta_ss.*rv_sd.*ra_dd.*rv_do).*ta_oo./(1-rv_dd.*ra_dd);
rtoa2           =   (ta_ss.*rv_sd+ta_sd.*rv_dd).*ta_do./(1-rv_dd.*ra_dd);
rtoa0           =   ra_so+ta_ss.*rv_so.*ta_oo;
R_TOA           =   T_g.*(rtoa0+rtoa1+rtoa2);
L_TOA           =   atm.La_.*R_TOA;

