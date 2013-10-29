%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains the properties of the evaporating compounds          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mixture of sub-cooled succinic acid adipic acid

% Molar gas constant (J/mol/K)
R = 8.314472;

% Number of species
nspec = 9;

% Molar masses of the compounds [kg/mol]
MW = [0.524 0.479 0.434    0.389   0.344   0.299 0.250 0.225 0.200];

% Volatilities of the surrogate compounds at T_ref (K). The volatilities can be
% given in either vapor pressures or saturation concentrations. 
T_ref = 298.15;
% % Vapor pressure in Pascals
% pstar = []; %
% Saturation concentration in ug/m3
%cstar = pstar.*MW.*1e9./R./T_ref 
cstar = [1e-5 1e-4 1e-3 1e-2 1e-1 1 1e1 1e2 1e3 ];
% Volatilities of the surrogate compounds at T_ref [Pa]
pstar = cstar*R.*T_ref./MW./1e9;

% Diffusion coefficients of the surrogate compounds in air at T_ref [m2/s]
% and the temperature-dependent factor
Dn = [5.0e-6 ] .* ones(1,nspec); 
mu1 = [1.75 ] .* ones(1,nspec);

% Densities of the surrogate compounds [kg/m3] 
rho = [1400 ] .* ones(1,nspec); % 

% Surface tensions of the surrogate compounds [N/m]
sigma1 = [0.05] .* ones(1,nspec); % 

% Vaporization enthalpies of the surrogate compounds at 298 K [J/mol]
%dHvap = 1e3.*[189.3134  169.5311  150.8227  133.1875  116.6248  101.1333 ]; % 
%dHvap = 1e3.*(85 - 11.*log10(cstar));
% dHvap = 1e3.*(131 - 11.*log10(cstar));
dHvap = dhvap_in .* ones(1,nspec) .* 1000;

% Mass and thermal accommodation coefficients of the surrogate species
alpha_m = [0.05 ] .* ones(1,nspec);
alpha_t = [1.0 ] .* ones(1,nspec);

% Initial mass fractions of the species in the aerosol
X_i = 1./nspec .* ones(1,nspec);



