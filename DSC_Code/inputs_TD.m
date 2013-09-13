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
%%%%%%%%%%%%%%%%%%  TD properties  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Length of the heating and cooling sections + total length of the TD [m]
l_heat = 0.55;
l_cool = 0.0;
l_tot = l_heat + l_cool;

% Residence time in the heating section [s]
t_res_heat = 14.0;
% Residence time in the cooling section [s]
t_res_cool = 0;
% Total residence time [s]
t_res_tot = t_res_heat + t_res_cool;

% Centerline velocity [m/s]
v = [];

if isempty(v)
    v = l_tot./t_res_tot;
elseif isempty(t_res_tot)
    t_res_tot = l_tot./v;
    t_res_heat = l_heat./v;
    t_res_cool = l_cool./v;
elseif isempty(l_tot)
    l_tot = v.*t_res_tot;
    l_heat = v.*t_res_heat;
    l_cool = v.*t_res_cool;
end



