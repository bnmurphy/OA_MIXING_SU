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
p_ref = 101325.0; 
% Volume [m3]
DSC.tot_vol = 80 .* 1.e-9;  %I think the DSC is about 80 mm3 
% DSC.tot_vol = DSC.liq_vol .* 10;  %Assume vapor has 10x the volume that the liquid does

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  DSC properties  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load Measurement Data (time, T, Q, Cp, Mass)
load([DSC_Data_Folder exp_name])

%Convert Measured Quantities
DSC.time = meas.time .* 60;  %min -> s
DSC.units.time = 'sec';
DSC.T = meas.T + 273.15; %deg C -> K
DSC.units.T = 'K';
DSC.Cp = meas.Cp .* 1000; %J/(g deg C) -> J/(kg K) 
DSC.units.Cp = 'J/(kg K)';
DSC.mass = meas.mass ./ 1e6; %mg -> kg
DSC.units.mass = 'kg';

%Confirm that Q is Defined the right way (positive for T increase)
switch DSC.T(end) > DSC.T(1)  
    case 1  %T increases throughout the experiment
        if sum(meas.Q) < 0, meas.Q = meas.Q .* -1; end
    case 0  %T decreases throughout the experiment
        if sum(meas.Q) > 0, meas.Q = meas.Q .* -1; end
end

DSC.Q = meas.Q .* 1000 .* DSC.mass;  %W/g -> W/kg -> W=J s-1
DSC.units.Q = 'J/s';
DSC.H2O_molfrac = meas.H2O_molfrac; %Mole fraction of Water Present

%Load initial Temperature
T_i = DSC.T(1);


DSC.mass = DSC.mass ./ 1;  %Reduce to particle amounts








