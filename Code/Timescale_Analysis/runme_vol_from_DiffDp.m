%Script for creating plots of timescale analysis
%  -comparing timescales of mixing within particles to evaporation
%  -X (dimensionless parameter) = tau_evap/tau_mix

clear variables

%% Setup Physical Parameters
pi = 3.1415;
R = 8.314; %J mol-1 K-1
T0 = 298;  %K
Dn_air = 5.0e-6;  %m2 s-1
mu1 = 1.75;  %Parameter for D_air calculation
MW = 0.1;  %kg mol-1
alpha_m = 1; %accomodation coefficient
rho1 = 1.4.*1000;  % Particle density (kg m-3)
sigma1 = 0.05; % Particle surface tension (N m-1)

Xm = 1.0;  %Assume mole fraction equal to one

%% Calculate Cutoff Satn Concentration as a Function of Diffusivity and Particle Size
%Mixing Diffusivity Range
Dmix_lo = 1.0e-20;  %m2 s-1
Dmix_hi = 1.0e-16;   %m2 s-1
Dmix = logspace(log10(Dmix_lo),log10(Dmix_hi),60);  %m2 s-1

%Size Range
Dp_lo = 1.0e-9; %m
Dp_hi = 10.0e-6; %m
Dp = logspace(-9, -5, 50); %m

%Mixing Equilibration Time (fxn of Dmix and Size)
tau_mix = zeros(length(Dmix),length(Dp));
for imix = 1:length(Dmix)
    for isize = 1:length(Dp)
        tau_mix(imix, isize) = calc_taumix_from_Dmix(Dp(isize),Dmix(imix)); %s
    end
end

%Evaporation Equilibration Time (=tau_mix), assume X = 1
tau_evap = tau_mix; %s

%Equilibrium Volatility associated with tau_evap
Dp = repmat(Dp,length(Dmix),1); %Array of particle radii
T = T0;
D_air = Dn_air.*(T/T0).^mu1;
Kn = calc_Knudsen(MW, D_air, Dp, T); %Calculate Knudsen Number
% Fuchs and Sutugin transition regime correction
Beta = calc_Beta(Kn, alpha_m);

csat = calc_csat_from_tauevap(tau_evap, T, Dp, D_air, Beta, Xm, MW, rho1, sigma1);

%Plot Diffusivity vs. Particle Size
x = Dmix .* 100^2; %cm2 s-1
y = Dp ;  %nm
z = log(csat)';
zticks = -5:6;
zlbl = -5:2:6;
lxlog = 1;
lylog = 1;
lsave = 1;
plot_contour(x, y, z, zticks, zlbl, 'Diffusivity (cm^2 s^{-1})', 'D_p (m)', ...
    'Volatility Required for Evaporation Timescale Equal to Mixing Timescale',...
    'none',[x(1),x(end)], [y(1),y(end)], 'Figs/', 'Vol_from_DiffDp', [], ...
    lxlog, lylog, lsave);









%% Temperature versus X for a series of volatilies




