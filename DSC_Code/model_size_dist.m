function mfr_dist = model_size_dist(X_i, c_aer_dist_int, alpha_in, TD_in, dhvap_in)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main program to model the dynamic evaporation in the TD     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all
% clear all

% Loading the inputs:
% Residence time, TD length and constants
inputs_TD
% Properties of the evaporating compounds
inputs_manish_modi
% Size distribution
inputs_size_dist_modi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%    Getting the composition of the aerosol and vapor mixture             %
%    at the initial temperature T_i. Aerosol is assumed to be internally  %
%    mixed.                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial mole fractions of the species in the aerosol
n_i_apu(1:nspec) = X_i./MW;
n_i_tot_apu = sum(X_i./MW);
Xm_i(1:nspec) = n_i_apu./n_i_tot_apu;

% Initial density of the aerosol
rhol_i = sum(X_i.*rho); % Mass-weighted average

% Initial surface tension of the aerosol
sigmal_i = sum(Xm_i.*sigma1); % Mole-weighted average

% Saturation pressures at initial temperature
psat_i(1:nspec) = pstar.*exp(dHvap.*(1./T_ref - 1./T_i)./R);

for i = 1:nspec
    %  Kelvin effect corresponding to the initial composition
    Ke_i(1:nbins+1,i) = exp(2.0.*MW(i).*sigmal_i./R./T_i./rho(i)./rp_i);
    % Mole fraction based equilibrium pressures
    %peq_i(1:nbins+1,i) = Xm_i(i).*psat_i(i).*Ke_i(1:nbins+1,i);
    % Mass fraction based equilibrium pressures
    peq_i(1:nbins+1,i) = X_i(i).*psat_i(i).*Ke_i(1:nbins+1,i);
end

% Initial partial pressures of the species, assuming aerosol initially in
% equilibrium with the aerosol corresponding to the peak size
pv_i(1:nspec) = peq_i(end,1:nspec);

% Initial particle mass
mp_i(1:nbins+1,1) = rhol_i.*4.0*pi.*rp_i.^3.0./3.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Calculation for each TD temperature                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:ntrials

Trial = k;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initializing the time-dependent variables at the TD temp for calculation %
%of evaporation in the heating section                                    %    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time
time0 = 0;
% Space
x_coord0 = l_heat.*time0./t_res_heat;

% if t_res_heat <17.0 && t_res_heat >13.0
% % Temperature profile for 16 s residence time in the beginning of the
% % reactor
%  para0 = [2.4509 -8.0918 5.0610 0.1405];
%  T_nd0 = para0(1).*x_coord0.^3 + para0(2).*x_coord0.^2 + para0(3).*x_coord0 + para0(4);
%  T_TD(1) = T_nd0.*(T_f(k) - T_i) + T_i;
% else
% % Flat temperature profile
T_TD(1) = T_f(k);
% end

% Equilibrium pressures at the TD temperature
psat(1:nspec) = pstar.*exp(dHvap.*(1./T_ref - 1./T_TD(1))./R);
csat(1:nspec) = MW.*psat./R./T_TD(1);


for i = 1:nspec
    % Kelvin effect corresponding to the initial composition
    Ke(1:nbins+1,i) = exp(2.0.*MW(i).*sigmal_i./R./T_TD(1)./rho(i)./rp_i);
    % Mole fraction based equilibrium pressures
    %peq(1:nbins+1,i) = Xm_i(i).*psat(i).*Ke(1:nbins+1,i);
    % Mass fraction based equilibrium pressures
    peq(1:nbins+1,i) = X_i(i).*psat(i).*Ke(1:nbins+1,i);
    % Partial pressure at the particle surface
    pv_a(1:nbins+1,i) = peq(1:nbins+1,i);
end

% Total particle mass concentration (kg/m3)
c_aer_tot = c_aer_tot_i.*T_i./T_TD(1);
% Particle mass concentrations (kg/m3)
c_aer_dist(1:nbins+1,1) = c_aer_dist_i.*T_i./T_TD(1);
% Particle number concentration (1/m3)
n_tot = n_tot_i.*T_i./T_TD(1);
% Particle number concentrations in each bin (1/m3)
n_dist(1:nbins+1,1) = n_dist_i.*T_i./T_TD(1);

% Diffusion coefficient (m2/s)
D(1:nspec) = Dn.*(T_TD(1)./T_ref).^mu1;

% Particle mass (kg) and size (m)
mp(1:nbins+1,1) = mp_i;
rp(1:nbins+1,1) = rp_i;

% Vapor mass concentrations (kg/m3)
cgas(1:nspec) = pv_i.*MW./R./T_TD(1);

% Total concentrations of the species 
ctot(1:nspec) = cgas + X_i.*c_aer_dist(end);
    
% Mean velocity of the gas molecules:
c_ave(1:nspec) = sqrt(8.*R.*T_TD(1)./MW./pi);
% Mean free path of the gas molecules:
lambda(1:nspec) = 3.*D./c_ave;

for i = 1:nspec
    % Masses of each species in the individual particles (kg)
    Pc(1:nbins+1,i) = X_i(i).*mp;
    % Knudsen number:
    Kn(1:nbins+1,i) = lambda(i)./rp_i;
    % Fuchs and Sutugin transition regime correction
    beta(1:nbins+1,i) = (1.0 + Kn(1:nbins+1,i))./...
        (1.0 + (4./alpha_m(i)./3 + 0.377).*Kn(1:nbins+1,i)...
        + 4.*Kn(1:nbins+1,i).^2./alpha_m(i)./3);
end

% Gas phase concentrations kg/m3 for the size distribution case and for the
% monodisperse case
Gc(1:2,1:nspec) = [cgas; cgas];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculating the time-dependent evaporation                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input(1:nbins+1,1:nspec) = Pc; % Masses of each species in each particle (kg)

input(nbins+2:nbins+3,1:nspec) = Gc; % Gas phase concentrations (kg/m3)


time = linspace(0,t_res_heat,100);
dt = mean(diff(time));

% Converting the input matrix to column vector for the odesolver
for i = 1:nspec
    input0(1+(i-1).*(nbins+3):i.*(nbins+3),1) = input(1:nbins+3,i); 
end


% Solving the mass fluxes 
options = odeset('RelTol',1E-6,'AbsTol',1E-19);    
[tout, output0] = ode45(@fluxes_size_dist, time, input0,  options, dt,nbins,nspec,l_heat,t_res_heat,...
    T_f(k),T_i,n_tot_i,n_dist_i,pstar,dHvap,T_ref,MW,sigma1,rho,Dn,mu1,p,alpha_m,alpha_t);

% % Converting back to our format
for i = 1:nspec
 output(1:length(time),1:nbins+3,i) = output0(1:length(time),1+(i-1).*(nbins+3):i.*(nbins+3));
end

% Temporary evolution of particle and gas phase masses (kg) and (kg/m3)
Pc_t_k(1:length(time),1:nbins+1,1:nspec) = output(1:length(time),1:nbins+1,1:nspec);
Gc_t_k(1:length(time),1:2,1:nspec) = output(1:length(time),nbins+2:nbins+3,1:nspec);

findex1 = find(Pc_t_k <= 0.0);
Pc_t_k(findex1) = 0.0; 
findex2 = find(Gc_t_k <= 0.0);
Gc_t_k(findex2) = 0.0; 

% Particle and gas phase masses at the end of the heating section (kg) and
% (kg/m3)
Pc_k_end(k,1:nbins+1,1:nspec) = Pc_t_k(end,1:nbins+1,1:nspec);
Gc_k_end(k,1:2,1:nspec) = Gc_t_k(end,1:2,1:nspec);

% Particle masses in the end (kg)
mp_end(k,1:nbins+1) = sum(Pc_k_end(k,1:nbins+1,1:nspec),3);

findex3 = find(mp_end <=0.0);
mp_end(findex3) = 0.0;

% Particle composition in the end
for i = 1:nspec
X_end(k,1:nbins+1,i) = Pc_k_end(k,1:nbins+1,i)./mp_end(k,1:nbins+1);
end

findex4 = find(X_end <= 0.0 | isnan(X_end));
X_end(findex4) = 0.0;

% Particle density (kg/m3)
for j = 1:nbins+1
         for i = 1:nspec
        apu(i) = X_end(k,j,i).*rho(i);
    end
    rhol_end(k,j) = sum(apu);
end

findex5 = find(rhol_end <= 0.0);
rhol_end(findex5) = 1000;

% Particle radii and diameters
rp_end(k,1:nbins+1) = (3.*mp_end(k,1:nbins+1)./4./pi./rhol_end(k,1:nbins+1)).^(1./3);
dp_end(k,1:nbins+1) = rp_end(k,1:nbins+1).*2;

% Final concentrations corrected back to T_i
% The number concentrations in each bin are the same as initially
c_aer_dist_end(k,1:nbins+1) = n_dist_i(1:nbins+1)'.*mp_end(k,1:nbins+1);
% Total mass
c_aer_tot_end(k) = sum(c_aer_dist_end(k,1:nbins));
% Mass fraction remaining
mfr_dist(k) = c_aer_tot_end(k)./c_aer_int_i;
mfr_mono(k) = c_aer_dist_end(k,nbins+1)./n_tot_i./mp_i(end);

org_conc_final = (c_aer_dist_end .* 1e9)';
org_conc_sum = sum(org_conc_final(1:nbins));

end
