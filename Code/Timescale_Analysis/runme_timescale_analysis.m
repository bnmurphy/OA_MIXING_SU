%Script for creating plots of timescale analysis
%  -comparing timescales of mixing within particles to evaporation
%  -X (dimensionless parameter) = tau_evap/tau_mix

clear variables

%% Setup Physical Parameters
pi = 3.1415;
R = 8.314; %J mol-1 K-1
T0 = 298;  %K
Dnair = 5.0e-6;  %m2 s-1
mu = 1.75;  %Parameter for D_air calculation
MW = 0.1;  %kg mol-1
alpha_m = 1; %accomodation coefficient
rho = 1.4.*1000;  % Particle density (kg m-3)
sigma = 0.05; % Particle surface tension (N m-1)


%%
%Mixing Diffusivity Range
Dmix_lo = 1.0e-20;  %m2 s-1
Dmix_hi = 2.0e-17;   %m2 s-1
Dmix = logspace(log(Dmix_lo),log(Dmix_hi),10);  %m2 s-1

%Size Range
Dp_lo = 1.0e-9; %m
Dp_hi = 10.0e-6; %m
Dp = [1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5]; %m

%Mixing Equilibration Time (fxn of Dmix and Size)
for imix = 1:length(Dmix)
    for isize = 1:length(Dp)
        tau_mix(imix, isize) = (Dp(isize)./2).^2 ./ pi.^2 ./ Dmix(imix); %s
    end
end

%Evaporation Equilibration Time (=tau_mix), assume X = 1
tau_evap = tau_mix; %s

%Equilibrium Volatility associated with tau_evap
rp = repmat(Dp./2,length(Dmix),1); %Array of particle radii
T = T0;
D_air = Dn_air.*(T/T0).^mu1(1:nspec);
c_ave = sqrt(8.*R.*T./MW./pi); % Mean velocity of the gas molecules
lambda = 3.*D_air./c_ave;      % Mean free path of the gas molecules
Kn = lambda./(Dp./2);          % Knudsen number
% Fuchs and Sutugin transition regime correction
Beta = (1.0 + Kn)./...
    (1.0 + (4./alpha_m./3 + 0.377).*Kn + 4.*Kn.^2./alpha_m./3);
peq_cutoff = R.*T0 ./ ...
    (Beta .* 4.*pi.* rp .* D_air .* tau_evap); %Pa

Xm = 1.0;  %Assume mole fraction equal to one
Ke = exp(2.0.*MW.*sigma./R./T_TD./rho./rp);  % Kelvin effect

% Saturation Pressure (Pa)
psat = peq ./ Xm ./ Ke;

%Saturation Concentration
csat = MW.*psat./R./T0;



%% Diffusivity versus cutoff volatility









%% Temperature versus X for a series of volatilies




