%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains the properties of the evaporating compounds          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mixture of sub-cooled succinic acid adipic acid

% Molar gas constant (J/mol/K)
R = 8.314472;

% Number of species
nspec = 10;

%Species in Mixture
spec_name = {'Malonic Acid',...
             'Succinic Acid',...
             'Glutaric Acid',...
             'Adipic Acid',...
             'Pimelic Acid',...
             'Suberic Acid',...
             'Azelaic Acid',...
             'Sebacic Acid',...
             'Dodecanoic Acid',...
             'Water'};

Xm_water = DSC.H2O_molfrac; %Mole fraction of water
iwater = 10;    %Index for water

%Specify Inputs for AIOMFAC - the activity coefficient model
%molec_group_flag is an array identifying the chemical identity of each
%funtional group present.
molec_group_flag = [2 43 0;...
                    2 43 0;...
                    2 43 0;...
                    2 43 0;...
                    2 43 0;...
                    2 43 0;...
                    2 43 0;...
                    2 43 0;...
                    1 2 43;...
                    17 0 0];
%molec_group_stoich is an array of the number of times each functional
%group (columns) appears in each molecule (row)
molec_group_stoich = [1 2 0;...
                      2 2 0;...
                      3 2 0;...
                      4 2 0;...
                      5 2 0;...
                      6 2 0;...
                      7 2 0;...
                      8 2 0;...
                      1 10 1;...
                      1 0 0];

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
% Dodecanoic acid pstar = 1.718e-9
pstar = [6.61e-9,...
         3.103e-9,...
         1.725e-9,...
         9.699e-10,...
         5.509e-10,...
         3.15e-10,...
         1.813e-10,...
         1.049e-10,...
         1.718e-9,...
         0.0413] .* 101325;
% Saturation concentration in ug/m3
cstar = pstar.*MW.*1e9./R./T_ref;

%Solubilitiy in Water  [mol solute / kg water]
Ksol = [13.45378,...
        0.49115,...
        3.25486,...
        0.15100,...
        0.56000,...
        0.01391,...
        0.00948,...
        0.00495,...
        0.00026,...
        1];

% Diffusion coefficients in air of the surrogate compounds in air at T_ref [m2/s]
% and the temperature-dependent factor
D_liq = [9.0e-10] .* ones(1,nspec-1);  %Diffusion of Organic Acids in Water
D_liq = [D_liq, 1e-3];  %Diffusion of Water in "sort of" water

Dn_air = [5.0e-6 ] .* ones(1,nspec-1); 
Dn_air = [Dn_air, 3.e-5]; %Add diffusion of water
mu1 = [1.75 ] .* ones(1,nspec);

% Liquid Densities of the surrogate compounds [kg/m3] 
rhol = [1.619,...
       1.56,...
       1.429,...
       1.36,...
       1.28,...
       1.272,...
       1.0287,...
       1.209,...
       0.88,...
       1.0] .* 1000;
rhos = rhol; %Let the solid density equal the liquid density for now.

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
     
%Fusion Enthalpies of the surrogate compounds at melting points [J/mol]
dHfus = [23.1,...
          34.0,...
          21.1,...
          33.7,...
          25.2,...
          41.8,...
          30.4,...
          47.0,...
          49.7,...
          6.01] .* 1000;
      
%Sublimation Enthalpy [J/mol]
dHsub = [132,...
        128,...
        134,...
        145,...
        153,...
        168,...
        178,...
        181,...
        169,...
        46.7].*1000;

%Fusion Entropy [J mol-1 K-1]
dSfus = [56.6,...
         74.7,...
         58,...
         80.4,...
         68.4,...
         101.2,...
         81.6,...
         116.4,...
         124.2,...
         22];

      
%Pure Component Melting Point at 1 atm [K]
Tmelt = [408,...
         455.2,...
         363.9,...
         419,...
         368.2,...
         413.2,...
         372.4,...
         403.9,...
         400.3,...
         273.15];

     
% Heat Capacity of Vapor, Liquid and Solid Phases
Cv_sld = [11,...
          14,...
          17,...
          20,...
          23,...
          26,...
          29,...
          32,...
          37] .*3.*R;  %J mol-1 K-1
Cv_sld = [Cv_sld, 37.94];  %Concatenate Solid Ice

Cp_liq = [linspace(220,650,nspec-1)]; %J mol-1 K-1
Cp_liq = [Cp_liq, 35.3]; %75.36  Concatenate Liquid Water

Cp_vap = 140 .* ones(1,nspec-1); %J mol-1 K-1
Cp_vap = [Cp_vap, 33.6]; %Concatenate Water Vapor
Cv_vap = Cp_vap; % - R;     %Constant Volume Heat Capacity

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
rhol_i = sum(X_i.*rhol); % Mass-weighted average

% Initial surface tension of the aerosol
sigmal_i = sum(Xm_i.*sigma1); % Mole-weighted average

%Calculate Initial DSC Properties Assuming STP
DSC.liq_vol = sum(X_i.*DSC.mass ./ rhol); %m3 if all mass is liquid

DSC.vap_vol = DSC.tot_vol - DSC.liq_vol; %Vapor volume m3
DSC.mass_conc = DSC.mass ./ DSC.tot_vol; %Total mass concentration (kg/m3)
DSC.L = DSC.vap_vol .^ (1/3);   %m - scaling for diffusion to the bulk
                                %    surface through the vapor
DSC.Area = DSC.L .^ 2; %m2
DSC.moles_air = p_ref .* DSC.vap_vol ./ R ./ T_ref; %[moles]
DSC.h2o = DSC.mass .* X_i(iwater); %kg of total water






