%%
% *************************************************************************
% *     functions: Run the SPART model for some exmaples                  *
% *     Authors:   Peiqi Yang (p.yang@utwente.nl)                         *
% *     Create:    10/Dec/2019                                          *
% *     Update:    10/MAY/2020                                            *
% *     Faculty of Geo-Information Science and Earth Observation (ITC)    *
% *     University of Twente, 7500 AE Enschede, The Netherlands           *
% *************************************************************************

% Soil-Plant-ATmosphere model (SPART) for top-of-canopy and top-of-atmosphere reflectance 
% Developed by Peiqi Yang (p.yang@utwente.nl),ITC, University of Twente
% Christiaan van der Tol and Wout Verhoef
% Coupling BSM, PROSAIL and SMAC to simulate TOA reflectance% 

% the example is in the folder [some_exameples]
clear 
clc
close all
restoredefaultpath
addpath('src_SPART','some_examples')          % add the path where the src is and will be removed in the end

% run Ex1_LSA_Cab_SPART.m
% run Ex2_LSA_LAI_SPART.m
run Ex3_LSA_AOT_SPART.m
% run Ex4_BRDF_SPART.m

restoredefaultpath