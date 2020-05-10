function [spectral] = define_bands

    % Define spectral regions for SPART v_1.20
    % All spectral regions are defined here as row vectors
    % originally by WV Jan. 2013
    % updated for SPART by PY Dec. 2019
    % 3 spectral regions for SPART
    
    reg1 =   400 :    1 :  2400;
    reg2 =  2500 :  100 : 15000;
    reg3 = 16000 : 1000 : 50000;
    
    spectral.wlS  = [reg1 reg2 reg3];
        
    % Other spectral (sub)regions
    
    spectral.wlP   = reg1;                            % PROSPECT data range
    spectral.wlE   = 400:1:750;                       % excitation in E-F matrix
    spectral.wlF   = 640:1:850;                       % chlorophyll fluorescence in E-F matrix
    spectral.wlO   = reg1;                            % optical part
    spectral.wlT   = [reg2 reg3];                     % thermal part   
    wlS            = spectral.wlS;                    % solar spectrum
    spectral.wlPAR = wlS(wlS>=400 & wlS<=700);        % PAR range
    
    nwlP            =   length(spectral.wlP);           % number of optical bands
    nwlT            =   length(spectral.wlT);           % number of thermal bands
    spectral.IwlP   =   1 : nwlP;                       % index of optical bands
    spectral.IwlT   =   nwlP+1 : nwlP+nwlT;             % index of thermal bands
    spectral.nwlP   =   nwlP;
    spectral.nwlT   =   nwlT;
end
