function [soilpar_i,leafbio_i,canopy_i,atm_i,angles_i] = select_parameter(soilpar,leafbio,canopy,atm,angles,k) 



soilpar_i         =  soilpar;
leafbio_i         =  leafbio;
canopy_i          =  canopy;
angles_i          =  angles;
atm_i             =  atm;

soilpar_i.B       = soilpar.B(k);
soilpar_i.lat     = soilpar.lat(k);
soilpar_i.lon     = soilpar.lon(k);
soilpar_i.SMp     = soilpar.SMp(k);

leafbio_i.Cab     = leafbio.Cab(k);           % chlorophyll content               [ug cm-2]
leafbio_i.Cdm     = leafbio.Cdm (k);          % dry matter content                [g cm-2]
leafbio_i.Cw      = leafbio.Cw(k);            % leaf water thickness equivalent   [cm]
leafbio_i.Cs      = leafbio.Cs(k);            % senescent material                [fraction]
leafbio_i.Cca     = leafbio.Cca(k);          % carotenoids                       [ug cm-2]
leafbio_i.Cant    = leafbio.Cant(k);         % anthocyanins                      [?]
leafbio_i.N       = leafbio.N(k);            % leaf structure parameter (affects the ratio of refl: transmittance) []

canopy_i.LAI      = canopy.LAI(k);
canopy_i.LIDFa    = canopy.LIDFa(k);
canopy_i.LIDFb    = canopy.LIDFb(k);
canopy_i.hot      = canopy.hot(k);           % ratio between leaf width: canopy height. This is the hot spot parameter. 
                                            % you were probably avoiding the hotspot,
                                            % so in that case the model is insensitive
                                            % to canopy.hot   
angles_i.tts      = angles.tts(k);           % solar zenith angle (deg)
angles_i.tto      = angles.tto(k);           % observation zenith angle (deg)
angles_i.psi      = angles.psi(k);           % absolute value (deg) of the difference between solar and zenith azimuth angle (arbitrary if angles.tto = 0)

atm_i.Pa          =  atm.Pa(k);              % Air pressure (hPa)
atm_i.aot550      =  atm.aot550(k);          % AOT at 550 nm
atm_i.uo3         =  atm.uo3 (k);            % Ozone content (cm)  0.3 cm= 300 Dobson Units
atm_i.uh2o        =  atm.uh2o(k);            % Water vapour (g/cm2)
atm_i.alt_m       =  atm.alt_m(k);           % altitude (m). In case, no air pressure input avabile
atm_i.Pa0         =  atm.Pa0 (k);            % sea level air pressure