function rad_conv   =    Spectral_convolution(wl_hi,radiation_spectra,sensor)
%% spectral convolution for a given spectral response function 
% input:
% wl_hi:        wavelength in high resoultion, to be convoluted
% radiation:    irradiance or radiance in high resolution, to be convoluted
% sensor:       [wl_srf][p_srf], wavelengts and conreponding contribution
%               wl_srf = number of bands contributed times number of bands
%               after convolution, wl_srf may contains nan value
%               p_srf, relative contribution. 
% Example:
%       wl_hi               =   400:500;
%       radiation_spectra   =   rand(length(wl_hi),1);
%       sensor.wl_srf       =   [400,401,402;450,451,452,453];
%       sensor.p_srf        =   [0.1,0.4,0.2;0.6,0.5,0.5,0.4];
%       Spectral_convolution(wl_hi,radiation_spectra,sensor)

wl_srf      =   sensor.wl_srf_smac;
p_srf       =   sensor.p_srf_smac;
[indx,~]    =   find_closest_all_elements(wl_srf,wl_hi);
rad         =   reshape(radiation_spectra(indx),size(wl_srf,1),size(wl_srf,2));
rad_conv    =   sum(rad.*p_srf,1)./sum(p_srf,1);

% NOTE, because p_srf has not been normalized, the normalization is needed.
% 
% wl_low      =   sensor.wl;
% indx_2      =   reshape(indx,size(wl_srf,1),size(wl_srf,2));


end 
function [closestIndex,closestValue] = find_closest_all_elements(V,N)
% V is the vectors you will go through 
% N is the target vector you want to find the minimal in. 
V       =   V(:);
N       =   N(:);
if size(N,2)~=1
    N   =   N';
end 
if size(V,2)~=1
    V   =   V';
end 
A                       =   repmat(N,[1 length(V)]);
[~,closestIndex]        =   min(abs(A-V'));
closestValue            =   N(closestIndex) ;

end 