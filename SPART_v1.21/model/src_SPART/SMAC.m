% SMAC model
function Opt_atm = SMAC(angles,atm)
%=============================================================================================
% library for atmospheric correction using SMAC method (Rahman and Dedieu, 1994)
% Translate by Peiqi Yang from C to Matlab 
% Improved by Peiqi Yang, ITC, University of Twente, from the original SMAC
%
% Peiqi Yang  (p.yang@utwente.nl)
% 26/Aug/2019
%

% use function Altitude2Pa


% Updates
% 21/Spet/2019 by Peiqi Yang
% - compute direct and diffuse transmittance for two paths
% - surface reflectance is not neccessary to be lambertian.
% - hot spot effects can be modelleded

% Variables:
% -input
% tts:          solar zenith angle
% tto:          observational zenith angle
% psi:          relative difference of solar and observational azimuth angle
% Pa:           Air pressure(can be estimated based on the elevation using Pa_altidude)
% aot550:       Aerosol optical thickness at 550 nm
% uO3:          Ozone content
% uh2o:         Water vapour content
% coef(a,n):    Fitting coefficients varying with spectral band
% TOC_R:        reflectance factors of surface, [rso,rsd,rdo,rdd], for
%               for lambertain surface, rso=rsd=rdo=rdd, which is the
%               orignial SMAC model assumes

% -important variable in the coding
% t+gas:        double path gaseous transmittance for each gas (to3,th2o,to2, tco2, tch4,tno2, tco)
% tg:           transtmittance for all gases
% Peq:          air press/ std air pressure, used to correction of air mass
% m, ms, mv:    two path air mass, air mass from incident to surface, and from
%               surface to TOA.
% taup:         aerosol optical depth in the spectral band, a function of
%               aot550. requires two fitting coefficients a0taup, a1taup

%=============================================================================================

if ~isempty(atm.Pa)||atm.Pa0~=-999
    atm.Pa = Altitude2Pa(atm.alt_m,atm.Pa0);
end
    
    
    tts         = angles.tts;
    tto         = angles.tto;
    psi         = angles.psi;
    Pa          = atm.Pa;
    taup550     = atm.aot550;
    uo3         = atm.uo3;
    uh2o        = atm.uh2o;
    coef        = atm.coef;
    
    ah2o    	= coef.ah2o;
    nh2o    	= coef.nh2o;
    ao3     	= coef.ao3;
    no3     	= coef.no3;
    ao2     	= coef.ao2;
    no2     	= coef.no2;
    po2     	= coef.po2;
    aco2        = coef.aco2;
    nco2        = coef.nco2;
    pco2        = coef.pco2;
    ach4        = coef.ach4;
    nch4        = coef.nch4;
    pch4        = coef.pch4;
    ano2        = coef.ano2;
    nno2        = coef.nno2;
    pno2        = coef.pno2;
    aco         = coef.aco;
    nco         = coef.nco;
    pco         = coef.pco;
    a0s         = coef.a0s;
    a1s         = coef.a1s;
    a2s         = coef.a2s;
    a3s         = coef.a3s;
    a0T         = coef.a0T;
    a1T         = coef.a1T;
    a2T         = coef.a2T;
    a3T     	= coef.a3T;
    taur    	= coef.taur;
    sr      	= coef.sr;
    a0taup  	= coef.a0taup;
    a1taup  	= coef.a1taup;
    wo      	= coef.wo;
    gc          = coef.gc;
    a0P     	= coef.a0P;
    a1P     	= coef.a1P;
    a2P     	= coef.a2P;
    a3P     	= coef.a3P;
    a4P     	= coef.a4P;
    Rest1   	= coef.Rest1;
    Rest2   	= coef.Rest2;
    Rest3   	= coef.Rest3;
    Rest4   	= coef.Rest4;
    Resr1   	= coef.Resr1;
    Resr2   	= coef.Resr2;
    Resr3   	= coef.Resr3;
    Resa1   	= coef.Resa1;
    Resa2   	= coef.Resa2;
    Resa3   	= coef.Resa3;
    Resa4   	= coef.Resa4;
    
    cdr         = pi/180;
    crd         = 180/pi;
    %--------------calculate the reflectance at TOA------------------------
    
    us = cos (tts * cdr);       % cos(tts)
    uv = cos (tto * cdr);       % cos(tto)
    Peq= Pa/1013.25;            % air press/ std air pressure
    
    %------------ 1) air mass  Eq(5) --------------------------------------
    m =  1/us + 1/uv;           %
    
    %----- 2) aerosol optical depth in the spectral band, taup Eq. (16)----
    taup = (a0taup) + (a1taup) * taup550 ;   %a0+a1*tau550
    
    
    %--3)gaseous transmissions (downward and upward paths)(Eq. 5 and 6)----
    %     to3 = 1;th2o= 1;to2 = 1; tco2= 1;tch4= 1;
    %   UH2O is the water vapour integrated content in g/cm2
    
    uo2 =  (Peq .^ (po2));      % Vertically integrated absorber amount
    uco2=  (Peq .^ (pco2));
    uch4=  (Peq .^ (pch4));
    uno2=  (Peq .^ (pno2));
    uco =  (Peq .^ (pco));
    
    % ao2(11:15) = ao2(11:15)*2;
    
    to3   = exp ( (ao3)  .* ( (uo3 *m)  .^ (no3)  ) );
    th2o  = exp ( (ah2o) .* ( (uh2o*m)  .^ (nh2o) ) );
    to2   = exp ( (ao2)  .* ( (uo2 *m)  .^ (no2)  ) );
    tco2  = exp ( (aco2) .* ( (uco2*m)  .^ (nco2) ) );
    tch4  = exp ( (ach4) .* ( (uch4*m)  .^ (nch4) ) );
    tno2  = exp ( (ano2) .* ( (uno2*m)  .^ (nno2) ) );
    tco   = exp ( (aco)  .* ( (uco *m)  .^ (nco) ) );
    
    tg      =   th2o .* to3 .* to2 .* tco2 .* tch4 .* tco .* tno2; %Eq. 6
    %%
    
    %-------5) spherical albedo of the atmosphere Eq. 7------------------------
    s = (a0s) * Peq +  (a3s) + (a1s)*taup550 + (a2s) * (taup550 .^ 2) ; % modification of Eq.8
    
    % -------6) Total scattering transmission Eq. 9----------------------------
    ttetas = (a0T) + (a1T)*taup550/us + ((a2T)*Peq + (a3T))/(1.+us);  % downward
    ttetav = (a0T) + (a1T)*taup550/uv + ((a2T)*Peq + (a3T))/(1.+uv);  % upward
    
    %----------7) scattering angle cosine  Eq.14-------------------------------
    cksi = - ( (us*uv) + (sqrt(1. - us*us) * sqrt (1. - uv*uv)*cos((psi) * cdr) ) );
    if cksi < -1        % it seems cksi is always>=-1;
        cksi=-1.0;
    end
    
    %--------------:  8) scattering angle in degree ---------------------------
    ksiD        =   crd*acos(cksi) ;
    %
    %--------------- 9) rayleigh atmospheric reflectance ----------------------
    ray_phase   =   0.7190443 * (1. + (cksi*cksi))  + 0.0412742;    %Eq. 13
    ray_ref     =   (taur*ray_phase ) / (4*us*uv);                  %Eq. 11
    ray_ref     =   ray_ref*Pa / 1013.25;                           %correction for pressure variation (not in the  paper)
    taurz       =   (taur)*Peq;                                     %Eq. 12
    
    % ------------ 10) aerosol atmospheric reflectance (3.4.2 )----------------
    aer_phase = a0P + a1P*ksiD + a2P*ksiD*ksiD +a3P*(ksiD .^ 3) + a4P * (ksiD .^ 4);    %extension of Eq. 17 aerosol phase function
    ak2 = (1. - wo).*(3. - wo.*3.*gc);                 % Eq. 15 k^2
    ak  = sqrt(ak2);                                % Eq. 15 k
    
    % ------------X Y Z Appendix--------------------------------------------------------
    e   = -3*us*us*wo ./  (4*(1 - ak2*us*us) );     %  E= -3*XMUS
    f   = -(1 - wo).*3.*gc.*us.*us.*wo ./ (4*(1 - ak2*us*us) );
    dp  = e / (3*us) + us*f;
    d   = e + f;
    b   = 2*ak ./ (3 - wo.*3.*gc);
    delta = exp( ak.*taup ).*(1 + b).^2 - exp(-ak.*taup).*(1 - b).^2 ;
    ww  = wo/4;
    ss  = us ./ (1 - ak2*us*us);
    q1  = 2 + 3*us + (1 - wo).*3.*gc*us*(1 + 2*us);
    q2  = 2 - 3*us - (1 - wo).*3.*gc*us*(1 - 2*us);
    q3  = q2.*exp( -taup/us );
    c1  =  ((ww.*ss) ./ delta) .* ( q1.*exp(ak.*taup).*(1 + b) + q3.*(1 - b) );
    c2  = -((ww.*ss) ./ delta) .* (q1.*exp(-ak.*taup).*(1 - b) + q3.*(1 + b) );
    cp1 =  c1.*ak ./ ( 3- wo.*3.*gc );
    cp2 = -c2.*ak ./ ( 3 - wo.*3.*gc );
    z   = d - wo.*3.*gc*uv.*dp + wo.*aer_phase/4;
    x   = c1 - wo.*3.*gc*uv.*cp1;
    y   = c2 - wo.*3.*gc*uv.*cp2;
    aa1 = uv ./ (1 + ak*uv);
    aa2 = uv ./ (1 - ak*uv);
    aa3 = us*uv / (us + uv) ;
    
    
    
    aer_ref1 = x.*aa1.* (1 - exp( -taup./aa1 ) )  ;              %Eq.13(1)
    aer_ref2 = y.*aa2.*( 1 - exp( -taup ./ aa2 )  ) ;            %Eq.13(2)
    aer_ref3 = z.*aa3.*( 1 - exp( -taup ./ aa3 )  ) ;               %Eq.13(3)
    % aer_ref3 = (z+aer_phase).*aa3.*( 1 - exp( -taup ./ aa3 )  ) ;    %Eq.13(3) in paper
    aer_ref = (aer_ref1+aer_ref2+aer_ref3) / ( us*uv );
    
    
    % ---------- 11) Residu Rayleigh(not in the paper)-------------------------
    Res_ray= Resr1 + Resr2 .* taur*ray_phase / (us*uv) + Resr3 .* ( (taur*ray_phase/(us*uv)) .^ 2) ;
    
    
    % ----------------:  12) Residu Aerosol-----------------------------------
    Res_aer= ( Resa1 + Resa2 .* ( taup * m *cksi ) + Resa3 .* ( (taup*m*cksi ) .^ 2) )+ Resa4 .* ( (taup*m*cksi).^ 3);
    
    %------------------ 13)  Term coupling molecule / aerosol--------------
    tautot=taup+taurz;
    
    Res_6s= ( Rest1+ Rest2 .* ( tautot * m *cksi )   + Rest3 .* ( (tautot*m*cksi) .^ 2) ) + Rest4 .* ( (tautot*m*cksi) .^ 3);
    
    % ----------------- 14) total atmospheric reflectance----------------------
    atm_ref = ray_ref - Res_ray + aer_ref - Res_aer + Res_6s;
    % --------------------15) TOA reflectance  --------------------------------
    % r_toa   = r_surf.*tg.*ttetas.*ttetav./(1-r_surf.*s) + (atm_ref .* tg) ;
    
    %% added by Peiqi Yang for non-lambertian surface
    tdir_tts          =     exp(-tautot/us);     % downward
    tdir_ttv          =     exp(-tautot/uv);     % upward
    tdif_tts          =     ttetas-tdir_tts;     % downward
    tdif_ttv          =     ttetav-tdir_ttv;     % upward
    
    Opt_atm.Ta_ss     =     tdir_tts;   % directional transmittance for direct incidence
    Opt_atm.Ta_sd     =     tdif_tts;   % hemispherical transmittance for direct incidence
    Opt_atm.Ta_oo     =     tdir_ttv;   % directional transmittance for direct incidence (in the viewing direction)
    Opt_atm.Ta_do     =     tdif_ttv;   % hemispherical transmittance for direct incidence(in the viewing direction)
    Opt_atm.Ta_s      =     ttetas;
    Opt_atm.Ta_o      =     ttetav;
    Opt_atm.Tg        =     tg;
    
    
    Opt_atm.Ra_dd     =     s;          % hemispherical atmospheric reflectance for diffuse light
    Opt_atm.Ra_so     =     atm_ref;    % directional atmospheric reflectance for direct incident
   
    function Pa  = Altitude2Pa(alt_m,Pa0,Ta0)
    %% Atmospheric pressure (in hpa) as a function of altitude (in meters)
    % Author: Peiqi Yang (p.yang@utwente.nl)
    %         28-Aug-2019
    
    % input: alt_m  =   altitude in meters
    %        Pa0    =   air pressure at sea level
    %        Ta0    =   air temperature at sea level
    % output: Pa= air pressure in hPa
    %
    % important constants:
    % normal temperature and pressure at sea level  Pa0 = 101325 (Pa)
    % temperature at sea level                      Ta0 = 15 degrees C or 288.15 degrees absolute
    % changes of temperature per 1km                dT  = 6.5 degrees per 1km
    % Earth-surface gravitational acceleration      g   = 9.80665 m/s2;
    % Molar mass of dry air                         M   = 0.02896968 kg/mol
    % Temperature lapse rate, = g/cp for dry air    L   = ~ 0.00976 K/m
    % Constant-pressure specific heat               Cp  = 1004.68506 J/(kg·K)
    % Universal gas constant                        R0  = 8.314462618 J/(mol·K)
    % Note: The gas constant is equivalent to the Boltzmann constant,
    %       but expressed in units of energy per temperature increment per mole,
    %       i.e. the pressure–volume product,
    %           rather than energy per temperature increment per particle.
    
    
    % [Barometric formula] (USED in here)
    % https://en.wikipedia.org/wiki/Barometric_formula
    % https://en.wikipedia.org/wiki/Atmospheric_pressure
    
    
    % Althernative: Standard Atmosphere model (USED in ORIGIANAL SMAC)
    % 1. Assumes that temperature is 15 degrees C at sea level (288.15 degrees absolute?
    % 2. Drops 6.5 degrees per 1000 meters of altitude, up to 11000 meters.
    % 3. It assumes that at sea level air pressure is 101325 Newtons per square meter,
    % 4. It asssumes air density is 1.225 kilograms per cubic meter.
    
    % Then
    % Pa/Pa0 = (Ta/Ta0)^2.558
    % Ta = Ta0-6.5/1000*alt_m;
    
    % reference: https://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html
    % https://www.quora.com/How-do-you-calculate-the-atmospheric-pressure-at-a-certain-altitude
    
    if nargin < 3
        Ta0 =   288.15;
    end
    if nargin < 2
        Pa0 =   1013.25;
    end
    
    g   =   9.80665;M   =   0.02896968; R0  =   8.314462618;
    Pa  =   Pa0.*exp(-(g*alt_m*M/(Ta0*R0)));                        %Pa_Altitude(1000) = 899.9491
    
    % dT  =   6.5; Cp  = 1004.68506;L   =   0.00976;
    % Pa  =   Pa0 * power( 1 - dT/1000 * alt_m / Ta0 , 5.25588);    %Pa_Altitude(1000) = 898.7456
    % 5.25588 can be changed to 5.31