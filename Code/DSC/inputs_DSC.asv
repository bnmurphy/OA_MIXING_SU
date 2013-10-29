%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This file contains the properties of the TD and other constants          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  Constants                           %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Molar gas constant [J/mol/K]
R = 8.314472;
% Boltzmann constant [J/K]
k = 1.3806503e-23; 
% Avogadro's number [1/mol]
Na = 6.0221415e23;
% Total pressure [Pa]
p = 101325.0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  DSC properties  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load Measurement Data (time, T, Q, Cp, Mass)
load('DSC_DATA/CappaMix_02.mat')

%Convert Measured Quantities
DSC.time = meas.time .* 60;  %min -> s
DSC.units.time = 'sec';
DSC.T = meas.T + 273.15; %deg C -> K
DSC.units.T = 'K';
DSC.Q = meas.Q .* 1000;  %W/g -> W/kg
DSC.units.Q = 'W/kg';
DSC.Cp = meas.Cp .* 1000; %J/(g deg C) -> J/(kg K) 
DSC.units.Cp = 'J/(kg K)';
DSC.mass = meas.mass ./ 1e6; %mg -> kg
DSC.units.mass = 'kg';

%Load initial Temperature
T_i = DSC.T(1);








