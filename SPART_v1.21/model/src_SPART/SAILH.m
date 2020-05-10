function [rad] = SAILH(soil,leafopt,canopy,angles)  
%% SAILH developed by W.V. 
%% edit by P.Y. (p.yang@utwente.nl)

% inputs
% soil.refl;            % [nwl]        soil reflectance spectra
% leafopt.refl;         % [nwl]        leaf/needle reflection  
% leafopt.tran;         % [nwl]        leaf/needle transmission
% canopy.LAI            % [0-8]        leaf area index
% canopy.LDIFa          % [-1,1]       leaf inclination distribution function
% canopy.LIDFb          % [-1,1]       leaf inclination distribution function 
% angles.tts;           %              solar zenith angle
% angles.tto;           %              observer zenith angle
% angles.psi;           %              relative azimuth angle

% output [Four reflectance factors] 
% rad.rsd     = rsd;
% rad.rdd     = rdd;
% rad.rdo     = rdo;
% rad.rso     = rso;


%% 
% pre-setting for the models
% Define canopy structure (SAIL assumpations)

canopy.nlayers      =   60;                             % number of layers
nl                  =   canopy.nlayers;                 % number of layers
canopy.x            =   (-1/nl : -1/nl : -1)';          % a column vector
canopy.xl           =   [0; canopy.x];                  % add top level
canopy.nlincl       =   13;                             % number of leaf inclination angle
canopy.nlazi        =   36;                             % number of leaf azimuth angle
canopy.litab        =   [5:10:75 81:2:89 ]';            % a column, never change the angles unless 'ladgen' is also adapted
canopy.lazitab      =   (5:10:355 );                    % a row     
canopy.lidf         =   leafangles(canopy.LIDFa,canopy.LIDFb);

% 0. Preparations

deg2rad = pi/180;

tts     = angles.tts;           % solar zenith angle
tto     = angles.tto;           % observer zenith angle
psi     = angles.psi;           % relative azimuth angle

nl      = canopy.nlayers;       % number of canopy layers (60)
litab   = canopy.litab;         % SAIL leaf inclibation angles
LAI     = canopy.LAI;           % leaf area index
lidf    = canopy.lidf;          % leaf inclination distribution function
x       = canopy.x;             % all levels except for the top
dx      = 1/nl;
iLAI    = LAI * dx;

rho     = leafopt.refl;         % [nwl]        leaf/needle reflection  
tau     = leafopt.tran;         % [nwl]        leaf/needle transmission
rs      = soil.refl;            % [nwl,nsoils] soil reflectance spectra
xl      = [0; x];               % [nl+1]       all levels + soil

%% 1.0 Geometric quantities

cos_tts     = cos(tts*deg2rad);             %           cos solar       angle   
tan_tto     = tan(tto*deg2rad);             %           tan observation angle

cos_tto     = cos(tto*deg2rad);             %           cos observation angle   
tan_tts     = tan(tts*deg2rad);             %           tan observation angle

psi         = abs(psi-360*round(psi/360));  %           (to ensure that volscatt is symmetric for psi=90 and psi=270)
dso         = sqrt(tan_tts.^2 + tan_tto.^2 - 2*tan_tts.*tan_tto.*cos(psi*deg2rad));

% 1.2 geometric factors associated with extinction and scattering
[chi_s,chi_o,frho,ftau]=volscat(tts,tto,psi,litab);   % volume scattering

cos_ttli    = cos(litab*deg2rad);           % [13]      cos leaf angles

ksli        = chi_s./cos_tts;               % [13]      p306{1} extinction coefficient in direction of sun        per leaf angle
koli        = chi_o./cos_tto;               % [13]      p307{1} extinction coefficient in direction of observer   per leaf angle

sobli       = frho*pi/(cos_tts*cos_tto);    % [13]      pag 309{1} area scattering coefficient fractions
sofli       = ftau*pi/(cos_tts*cos_tto);    % [13]      pag 309{1}
bfli        = cos_ttli.^2;                  % [13]

%integration over angles (using a vector inproduct) -> scalars
k           = ksli'*lidf;                   %           pag 306{1}    extinction coefficient in direction of sun.
K           = koli'*lidf;                   %           pag 307{1}    extinction coefficient in direction of observer
bf          = bfli'*lidf;                   % 
sob         = sobli'*lidf;                  %           weight of specular2directional back    scatter coefficient
sof         = sofli'*lidf;                  %           weight of specular2directional forward scatter coefficient
% 1.3 geometric factors to be used later with rho and tau, f1 f2 of pag 304:
% these variables are scalars 
sdb         = 0.5*(k+bf);                   % fs*f1
sdf         = 0.5*(k-bf);                   % fs*f2     weight of specular2diffuse     foward  scatter coefficient 
ddb         = 0.5*(1+bf);                   % f1^2+f2^2 weight of diffuse2diffuse      back    scatter coefficient 
ddf         = 0.5*(1-bf);                   % 2*f1*f2   weight of diffuse2diffuse      forward scatter coefficient 
dob         = 0.5*(K+bf);                   % fo*f1     weight of diffuse2directional  back    scatter coefficient 
dof         = 0.5*(K-bf);                   % fo*f2     weight of diffuse2directional  forward scatter coefficient 

% 1.5 probabilities Ps, Po, Pso
Ps          =   exp(k*xl*LAI);                                              % [nl+1]  p154{1} probability of viewing a leaf in solar dir
Po          =   exp(K*xl*LAI);                                              % [nl+1]  p154{1} probability of viewing a leaf in observation dir

Ps(1:nl)    =   Ps(1:nl) *(1-exp(-k*LAI*dx))/(k*LAI*dx);                    % Correct Ps/Po for finite dx
Po(1:nl)    =   Po(1:nl) *(1-exp(-K*LAI*dx))/(K*LAI*dx);                    % Correct Ps/Po for finite dx

q           =   canopy.hot;
Pso         =   zeros(size(Po));
for j=1:length(xl)
    Pso(j,:)=   quad(@(y)Psofunction(K,k,LAI,q,dso,y),xl(j)-dx,xl(j))/dx;
end
% keyboard
% q           =   canopy.hot;
% alf         =   dso/q*2./(k+K);
% 
% dPso        =   [1,  1-(K+k-sqrt(K.*k).*(exp(-alf.*(xl(1:60)'+0.9833)*iLAI)))*iLAI];  Pso= cumprod(dPso);
% Pso(1:nl)   =   0.5*(Pso(1:nl)+Pso(2:nl+1));        % a simple way to correct Pso for the finite dx 

Pso(Pso>Po)= min([Po(Pso>Po),Ps(Pso>Po)],[],2);    %takes care of rounding error
Pso(Pso>Ps)= min([Po(Pso>Ps),Ps(Pso>Ps)],[],2);    %takes care of rounding error

% the following are vectors with lenght nwl
sigb        = ddb*rho + ddf*tau;            % [nwl]     sigmab, p305{1} diffuse     backscatter scattering coefficient for diffuse  incidence 
sigf        = ddf*rho + ddb*tau;            % [nwl]     sigmaf, p305{1} diffuse     forward     scattering coefficient for forward  incidence 
sb          = sdb*rho + sdf*tau;            % [nwl]     sb,     p305{1} diffuse     backscatter scattering coefficient for specular incidence 
sf          = sdf*rho + sdb*tau;            % [nwl]     sf,     p305{1} diffuse     forward     scattering coefficient for specular incidence 
vb          = dob*rho + dof*tau;            % [nwl]     vb,     p305{1} directional backscatter scattering coefficient for diffuse  incidence 
vf          = dof*rho + dob*tau;            % [nwl]     vf,     p305{1} directional forward     scattering coefficient for diffuse  incidence 
w           = sob*rho + sof*tau;            % [nwl]     w,      p309{1} bidirectional scattering coefficent (directional-directional)         
a           = 1-sigf;                       % [nwl]     attenuation
m           = sqrt(a.^2-sigb.^2);           % [nwl]
rinf        = (a-m)./sigb;                  % [nwl]
rinf2       = rinf.*rinf;                   % [nwl]

% direct solar radiation
J1k        = calcJ1(-1, m,k,LAI);          % [nwl]
J2k        = calcJ2( 0, m,k,LAI);          % [nwl]
J1K        = calcJ1(-1, m,K,LAI);          % [nwl]   % added for calculation of rdo
J2K        = calcJ2( 0, m,K,LAI);          % [nwl]   % added for calculation of rdo

e1          = exp(-m*LAI);                  % [nwl]
e2          = e1.^2;                        % [nwl]
re          = rinf.*e1;                     % [nwl]

denom       = 1-rinf2.*e2;                  % [nwl]

s1          = sf+rinf.*sb;
s2          = sf.*rinf+sb;
v1          = vf+rinf.*vb;
v2          = vf.*rinf+vb;

Pss         = s1.*J1k;          % [nwl]
Qss         = s2.*J2k;          % [nwl]

Poo         = v1.*J1K;          % (nwl)   % added for calculation of rdo
Qoo         = v2.*J2K;          % [nwl]   % added for calculation of rdo

tau_ss      = exp(-k*LAI);                  % [1]
tau_oo      = exp(-K*LAI);                  % [1]

Z           = (1 - tau_ss * tau_oo)/(K + k);  % needed for analytic rso

tau_dd      = (1-rinf2).*e1 ./denom;        % [nwl]
rho_dd      = rinf.*(1-e2)  ./denom;        % [nwl]
tau_sd      = (Pss-re.*Qss) ./denom;        % [nwl]
tau_do      = (Poo-re.*Qoo) ./denom;        % [nwl]
rho_sd      = (Qss-re.*Pss) ./denom;        % [nwl]
rho_do      = (Qoo-re.*Poo) ./denom;        % (nwl)

T1          = v2.*s1.*(Z-J1k*tau_oo)./(K+m)+v1.*s2.*(Z-J1K*tau_ss)./(k+m);
T2          = -(Qoo.*rho_sd+Poo.*tau_sd).*rinf;
rho_sod     = (T1+T2)./(1-rinf2);

rho_sos     = w * sum(Pso(1:nl))*iLAI;
rho_so      = rho_sod + rho_sos;

Pso2w       = Pso(nl+1);

% SAIL analytical reflectances

denom       =   1-rs.*rho_dd;

rso         =   rho_so + rs * Pso2w                                        ...
              + ((tau_sd+tau_ss*rs.*rho_dd)*tau_oo+(tau_sd+tau_ss).*tau_do) ...
               .*rs./denom; 
           
rdo         =   rho_do + (tau_oo + tau_do).*rs.*tau_dd./denom;        

rsd         =   rho_sd + (tau_ss + tau_sd).*rs.*tau_dd./denom;
rdd         =   rho_dd + tau_dd.*rs.*tau_dd./denom;



rad.rsd     = rsd;
rad.rdd     = rdd;
rad.rdo     = rdo;
rad.rso     = rso;
return

%% APPENDIX I functions J1 and J2 (introduced for numerically stable solutions)

function J1 = calcJ1(x,m,k,LAI)

    % Modified by W. Verhoef on 14 Jan 2016 to repair a malfunctioning of
    % RTMo_lite; should be inserted in RTMo as well

    J1 = zeros(length(m),1);

    sing = abs((m-k)*LAI) < 1e-6;   % near singularity logical array
    
    CS = find(sing);                % indices singular case
    CN = find(~sing);               % indices normal case
        
    J1(CN) = (exp(m(CN)*LAI*x)-exp(k*LAI*x))./(k-m(CN));
    
    J1(CS) = -.5*(exp(m(CS)*LAI*x)+exp(k*LAI*x))*LAI.*x.*(1-1/12*(k-m(CS)).^2*LAI^2.*x.^2);

return

function J2 = calcJ2(x,m,k,LAI)
    J2 = (exp(k*LAI*x)-exp(-k*LAI)*exp(-m*LAI*(1+x)))./(k+m);
return;

%% APPENDIX II function volscat

function [chi_s,chi_o,frho,ftau]    =   volscat(tts,tto,psi,ttli)

%Volscatt version 2.
%created by W. Verhoef
%edited by Joris Timmermans to matlab nomenclature.
% date: 11 February 2008
%tts    [1]         Sun            zenith angle in degrees
%tto    [1]         Observation    zenith angle in degrees
%psi    [1]         Difference of  azimuth angle between solar and viewing position
%ttli   [ttli]      leaf inclination array

deg2rad = pi/180;
nli     = length(ttli);

psi_rad         = psi*deg2rad*ones(nli,1);

cos_psi         = cos(psi*deg2rad);                 %   cosine of relative azimuth angle

cos_ttli        = cos(ttli*deg2rad);                %   cosine of normal of upperside of leaf
sin_ttli        = sin(ttli*deg2rad);                %   sine   of normal of upperside of leaf

cos_tts         = cos(tts*deg2rad);                 %   cosine of sun zenith angle
sin_tts         = sin(tts*deg2rad);                 %   sine   of sun zenith angle

cos_tto         = cos(tto*deg2rad);                 %   cosine of observer zenith angle
sin_tto         = sin(tto*deg2rad);                 %   sine   of observer zenith angle

Cs              = cos_ttli*cos_tts;                 %   p305{1}
Ss              = sin_ttli*sin_tts;                 %   p305{1}

Co              = cos_ttli*cos_tto;                 %   p305{1}
So              = sin_ttli*sin_tto;                 %   p305{1}

As              = max([Ss,Cs],[],2);
Ao              = max([So,Co],[],2);

bts             = acos(-Cs./As);                    %   p305{1}
bto             = acos(-Co./Ao);                    %   p305{2}

chi_o           = 2/pi*((bto-pi/2).*Co + sin(bto).*So);
chi_s           = 2/pi*((bts-pi/2).*Cs + sin(bts).*Ss);

delta1          = abs(bts-bto);                     %   p308{1}
delta2          = pi-abs(bts + bto - pi);           %   p308{1}

Tot             = psi_rad + delta1 + delta2;        %   pag 130{1}

bt1             = min([psi_rad,delta1],[],2);
bt3             = max([psi_rad,delta2],[],2);
bt2             = Tot - bt1 - bt3;

T1              = 2.*Cs.*Co + Ss.*So.*cos_psi;
T2              = sin(bt2).*(2*As.*Ao + Ss.*So.*cos(bt1).*cos(bt3));

Jmin            = (   bt2).*T1 - T2;
Jplus           = (pi-bt2).*T1 + T2;

frho            =  Jplus/(2*pi^2);
ftau            = -Jmin /(2*pi^2);

% pag.309 wl-> pag 135{1}
frho            = max([zeros(nli,1),frho],[],2);
ftau            = max([zeros(nli,1),ftau],[],2);
return

%% APPENDIX IV function Pso

function pso    =   Psofunction(K,k,LAI,q,dso,xl)
if dso~=0
    alf         =   (dso/q) *2/(k+K);
    
    pso         =   exp((K+k)*LAI*xl + sqrt(K*k)*LAI/(alf  )*(1-exp(xl*(alf  ))));% [nl+1]  factor for correlation of Ps and Po
else
    pso         =   exp((K+k)*LAI*xl - sqrt(K*k)*LAI*xl);% [nl+1]  factor for correlation of Ps and Po
end

%% APPENDIX IV function LIDF
function [lidf]=  leafangles(a,b)                                     
% Subroutine FluorSail_dladgen
% Version 2.3 
% For more information look to page 128 of "theory of radiative transfer models applied in optical remote sensing of
% vegetation canopies"
%
% FluorSail for Matlab
% FluorSail is created by Wout Verhoef, 
% National Aerospace Laboratory (NLR)
% Present e-mail: w.verhoef@utwente.nl
%
% This code was created by Joris Timmermans, 
% International institute for Geo-Information Science and Earth Observation. (ITC)
% Email: j.timmermans@utwente.nl
%
%% main function
F           =   zeros(1,13);
for i=1:8                                                               
    theta   =   i*10;                  %                theta_l =  10:80
    F(i)    =   dcum(a,b,theta);      %                 F(theta)
end

for i=9:12                                                              
    theta   =   80 + (i-8)*2;                         % theta_l = 82:88
    F(i)    =   dcum(a,b,theta);                     %  F(theta)
end

for i=13:13                                          %  theta_l = 90:90
    F(i) =   1;                                      %  F(theta)
end

lidf        =   zeros(13,1);
for i=13:-1:2                                                           
    lidf(i) =   F(i) -   F(i-1);                     %  lidf   =   dF/dtheta;
end
lidf(1) =   F(1);                                    %  Boundary condition

%% SubRoutines
function [F]   =  dcum(a,b,theta)
rd  =   pi/180;                                     %   Geometrical constant
if a>1 
    F    =   1 - cos(theta*rd);
else
    eps     =   1e-8;
    delx    =   1;
    
    x       =   2*rd *theta;
    theta2  =   x;
                                                                        %    
    while (delx > eps)
        y   =   a*sin(x) + 0.5*b*sin(2*x);
        dx  =   0.5*(y - x + theta2);
        x   =   x + dx;
        delx=   abs(dx);
    end
    F    =   (2*y + theta2)/pi;                     %   Cumulative leaf inclination density function
    %pag 139 thesis says: F = 2*(y+p)/pi. 
    %Since theta2=theta*2 (in rad), this makes F=(2*y + 2*theta)/pi    
end