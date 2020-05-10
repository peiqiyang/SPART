%%
% *************************************************************************
% *     functions: read SMAC coefficents into a structure two-way         *
% *     Authors:   Peiqi Yang (p.yang@utwente.nl)                         *
% *     Date:      10/March/2020                                          *
% *     Faculty of Geo-Information Science and Earth Observation (ITC)    *
% *     University of Twente, 7500 AE Enschede, The Netherlands           *
% *************************************************************************

% read SMAC coefficents in .dat files into a structure
% input: filename
% output: coef
% filename = 'SMAC_Sentinel3\coef_SENTINEL3A_OLCI_V32_CONTINENTAL.dat';

function coef = read_coeffs(filename)

datafile    =   importdata(filename);
if isstruct(datafile)
    coefs       =   datafile.data;
else
    coefs = datafile;
end
% nwl         =   size(coefs,2);
ncoefs        = size(coefs,1);
if ncoefs ~=   49
    coefs   =   coefs';
end
%%  coefficients for gaseous transmission
% H2O
coef.ah2o   =   coefs(1,:);
coef.nh2o   =   coefs(2,:);
% O3
coef.ao3   =   coefs(3,:);
coef.no3   =   coefs(4,:);
% O2
coef.ao2   =   coefs(5,:);
coef.no2   =   coefs(6,:);
coef.po2   =   coefs(7,:);
% CO2
coef.aco2   =   coefs(8,:);
coef.nco2   =   coefs(9,:);
coef.pco2   =   coefs(10,:);
% NH4
coef.ach4   =   coefs(11,:);
coef.nch4   =   coefs(12,:);
coef.pch4   =   coefs(13,:);

% NO2
coef.ano2   =   coefs(14,:);
coef.nno2   =   coefs(15,:);
coef.pno2   =   coefs(16,:);

% CO
coef.aco   =   coefs(17,:);
coef.nco   =   coefs(18,:);
coef.pco   =   coefs(19,:);

% rayleigh scattering (section 3.3)
coef.a0s   =   coefs(20,:);
coef.a1s   =   coefs(21,:);
coef.a2s   =   coefs(22,:);
coef.a3s   =   coefs(23,:);
% aerosol scattering (section 3.3)
coef.a0T   =   coefs(24,:);
coef.a1T   =   coefs(25,:);
coef.a2T   =   coefs(26,:);
coef.a3T   =   coefs(27,:);
%
coef.taur   =   coefs(28,:);
coef.sr     =   coefs(29,:);

coef.a0taup  =  coefs(30,:);
coef.a1taup  =  coefs(31,:);


coef.wo  =  coefs(32,:);
coef.gc  =  coefs(33,:);

% coefficients for aersol phase function
coef.a0P  =  coefs(34,:);
coef.a1P  =  coefs(35,:);
coef.a2P  =  coefs(36,:);

% coefficients for the reseduals
coef.a3P  =  coefs(37,:);
coef.a4P  =  coefs(38,:);

coef.Rest1  =  coefs(39,:);
coef.Rest2  =  coefs(40,:);

coef.Rest3  =  coefs(41,:);
coef.Rest4  =  coefs(42,:);


coef.Resr1  =  coefs(43,:);
coef.Resr2  =  coefs(44,:);
coef.Resr3  =  coefs(45,:);

coef.Resa1  =  coefs(46,:);
coef.Resa2  =  coefs(47,:);

coef.Resa3  =  coefs(48,:);
coef.Resa4  =  coefs(49,:);





