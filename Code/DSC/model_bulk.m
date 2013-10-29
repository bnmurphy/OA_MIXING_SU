function mfr_dist = model_size_dist(X_i, c_aer_dist_int, alpha_in, TD_in, dhvap_in)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main program to model the dynamic evaporation in the        %
%   Differential Scanning Calorimeter. An energy profile and mixture      %
%   properties are input to the model and temperature and concentration   %
%   profiles are output
%
% Although the DSC experiment focuses on bulk fluid, try to keep particle
% functionality throughout the code, so that it can be extended to look at
% particle issues afterward
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all
% clear all

% Loading the inputs:
% DSC properties and Experimental Data
inputs_DSC
% Properties of the evaporating compounds
inputs_chem

% Size distribution
%  (Comment out now because the DSC system treats one flat surface)
% inputs_aero

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%    Getting the composition of the aerosol and vapor mixture             %
%    at the initial temperature T_i. Aerosol is assumed to be internally  %
%    mixed.                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Saturation pressures at initial temperature [Pa]
psat_i(1:nspec) = pstar.*exp(dHvap.*(1./T_ref - 1./T_i)./R);

for i = 1:nspec
    % Calculate Kelvin effect
    %Ke_i(1:nbins+1,i) = exp(2.0.*MW(i).*sigmal_i./R./T_i./rho(i)./rp_i);
    %  No Kelvin effect for DSC -> flat surface
    Ke_i(i) = 1.0;
    
    % Mole fraction based equilibrium pressures
    %peq_i(1:nbins+1,i) = Xm_i(i).*psat_i(i).*Ke_i(1:nbins+1,i);
    
    % Mass fraction based equilibrium pressures
    peq_i(i) = X_i(i).*psat_i(i).*Ke_i(i);
end

% Initial partial pressures of the species, assuming aerosol initially in
% equilibrium with the aerosol corresponding to the peak size
pv_i(1:nspec) = peq_i(end,1:nspec);

% Initial particle mass
mp_i(1:nspec) = X_i .* DSC.mass; %kg


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initializing the time-dependent variables                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Equilibrium pressures at the DSC temperature
psat = psat_i;
csat(1:nspec) = MW.*psat./R./T_i;

%Equilibrium Pressure initially
peq = peq_i;
pv_a = pv_i;  %Initial Partial pressure at equilibrium

% Vapor mass concentrations (kg/m3), initial
cgas(1:nspec) = pv_i.*MW./R./T_i;

% Diffusion coefficient (m2/s)
D(1:nspec) = Dn.*(T_i./T_ref).^mu1;

% Bulk mass concentration (kg/m3)
Bc = mp_i;


% Total concentrations of the species 
% ctot(1:nspec) = cgas + X_i.*c_aer_dist(end);
ctot(1:nspec) = cgas + X_i.*mp(1,1);
    
% Mean velocity of the gas molecules:
c_ave(1:nspec) = sqrt(8.*R.*T_i./MW./pi);
% Mean free path of the gas molecules:
lambda(1:nspec) = 3.*D./c_ave;


% Gas phase concentrations kg/m3 for the size distribution case and for the
% monodisperse case
Gc(1:nspec) = cgas;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculating the time-dependent evaporation                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input = [T_i, Bc, Gc]'; % T, Temperature, K  
                        % Bc, Masses of each species bulk (kg)
                        % Gc, Gas phase concentrations (kg/m3)
time = DSC.time;
dt = mean(diff(time));

% Solving the mass fluxes 
options = odeset('RelTol',1E-6,'AbsTol',1E-19);    
[tout, output] = ode45(@fluxes_bulk, time, input0,  options, dt,nspec,...
    pstar,dHvap,T_ref,MW,sigma1,rho,Dn,mu1,p,alpha_m);

T_out = output(:,1);
Bc_out = output(:,2:nspec+1);
Gc_out = output(:,nspec+2, 2*nspec+1);




end
