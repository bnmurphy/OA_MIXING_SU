%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains the properties of the evaporating compounds          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mixture of sub-cooled succinic acid adipic acid

% Molar gas constant (J/mol/K)
R = 8.314472;

% Number of species
nspec = 10;

%Species in Mixture
% Malonic Acid
% Succinic Acid
% Glutaric Acid
% Adipic Acid
% Pimelic Acid
% Suberic Acid
% Azelaic Acid
% Sebacic Acid
% Dodecanoic Acid
% Water
Xm_water = DSC.H2O_molfrac; %Mole fraction of water
iwater = 10;    %Index for water

% Molar masses of the compounds [kg/mol]
MW = [104.06,...
      118.09,...
      132.11,...
      146.14,...
      160,...
      174,...
      188,...
      202,...
      230,...
      18] ./ 1000;

% Volatilities of the surrogate compounds at T_ref (K). The volatilities can be
% given in either vapor pressures or saturation concentrations. 
T_ref = 298.15;
% % Vapor pressure in Pascals
pstar = [6.61e-9,...
         3.103e-9,...
         1.725e-9,...
         9.699e-10,...
         5.509e-10,...
         3.15e-10,...
         1.813e-10,...
         1.049e-10,...
         1.718e-9,...
         0.023024] .* 101325;
% Saturation concentration in ug/m3
cstar = pstar.*MW.*1e9./R./T_ref;

% Diffusion coefficients in air of the surrogate compounds in air at T_ref [m2/s]
% and the temperature-dependent factor
Dn = [5.0e-6 ] .* ones(1,nspec); 
mu1 = [1.75 ] .* ones(1,nspec);

% Densities of the surrogate compounds [kg/m3] 
rho = [1.619,...
       1.56,...
       1.429,...
       1.36,...
       1.28,...
       1.272,...
       1.0287,...
       1.209,...
       0.88,...
       1] .* 1000;

% Surface tensions of the surrogate compounds [N/m]
sigma1 = [0.05] .* ones(1,nspec); % 

% Vaporization enthalpies of the surrogate compounds at 298 K [J/mol]
dHvap = [141.9,...
         88.5,...
         141.0,...
         111.0,...
         124.0,...
         130.0,...
         146.0,...
         140.0,...
         119.0,...
         40.14] .* 1000;

% Heat Capacity of Vapor and Liquid Phases
Cp_liq = linspace(330,450,nspec-1); %J mol-1 K-1
Cp_liq = [Cp_liq, 75.3]; %Concatenate Liquid Water
Cp_vap = 140 .* ones(1,nspec-1); %J mol-1 K-1
Cp_vap = [Cp_vap, 33.6]; %Concatenate Liquid Water

% Mass and thermal accommodation coefficients of the surrogate species
alpha_m = [0.05 ] .* ones(1,nspec);
alpha_t = [1.0 ] .* ones(1,nspec);


% Initial mole fractions of the species in the condensed phase
Xm_i(iwater) = Xm_water;
Xm_i(1:nspec-1) = ones(1,nspec-1) ./(nspec-1) .* (1-Xm_water); %Assume equi-molar

% Initial mass fractions of the species in the condensed phase
m_i_apu(1:nspec) = Xm_i.*MW;
m_i_tot_apu = sum( m_i_apu );
X_i(1:nspec) = m_i_apu./m_i_tot_apu;

% Initial density of the aerosol
rhol_i = sum(X_i.*rho); % Mass-weighted average

% Initial surface tension of the aerosol
sigmal_i = sum(Xm_i.*sigma1); % Mole-weighted average

%Calculate Initial Mass concentration
DSC.liq_vol = DSC.mass ./ rhol_i; %m3
DSC.tot_vol = DSC.liq_vol .* 10;  %Assume vapor has 10x the volume that the liquid does
DSC.vap_vol = DSC.tot_vol - DSC.liq_vol; %Vapor volume
DSC.mass_conc = DSC.mass ./ DSC.tot_vol; %Total mass concentration (kg/m3)
DSC.area = DSC.tot_vol .^ 2/3;  %m2




















