%Script for creating plots of timescale analysis
%  -comparing timescales of mixing within particles to evaporation
%  -X (dimensionless parameter) = tau_evap/tau_mix

clear variables

%% Setup Physical Parameters
pi = 3.1415;
R = 8.314; %J mol-1 K-1
k = 1.38e-23; %m2 kg s-2 K-1
T0 = 298;  %K
Dn_air = 5.0e-6;  %m2 s-1
mu1 = 1.75;  %Parameter for D_air calculation
rh = 2.12e-10; %Hydrodynamic Radius (m)
MW = 0.1;  %kg mol-1
alpha_m = 1; %accomodation coefficient
rho1 = 1.4.*1000;  % Particle density (kg m-3)
sigma1 = 0.05; % Particle surface tension (N m-1)
dHvap = 100*1000; %J mol-1

Xm = 1.0;  %Assume mole fraction equal to one


%% Set Independent Variable Ranges
%Temperature Range
T_lo = 273.15; %K
T_hi = 90+273.15; %K
% T_range = linspace(T_lo,T_hi,1000);
T_range = [0 25 90]+273.15;

%Water Composition Range
H2O_lo = 0.0;
H2O_hi = 0.95;
H2O_range = linspace(H2O_lo, H2O_hi, 1200);
% H2O_range = [0.0 0.25 0.50 0.75 0.9 0.95];

%Size Range
Dp_lo = 1.0e-9; %m
Dp_hi = 10.0e-6; %m
Dp = logspace(-9, -5, 1000); %m

%% Find Diffusivity - interpolate H2O and T
%Load Viscosity Measurements
load('expt_visc.mat')


%% Find Diffusivity - interpolate H2O and T
H2Omod = repmat(H2O_range,length(T_range),1);
Tmod = repmat(T_range',1,length(H2O_range));
visc = 10.^ interp2(expt.T'+273.15, expt.water,...
    log10(expt.visc.real'), Tmod, H2Omod,'spline');
Dmix = k .* Tmod ./ 6 ./ pi ./ visc ./ rh; %m2 s-1


%% Loop through sizes
cstar = zeros(length(T_range),length(H2O_range),length(Dp));
for isize = 1:length(Dp)
    
    %Mixing Equilibration Time (fxn of Dmix and Size)
    tau_mix = calc_taumix_from_Dmix(Dp(isize),Dmix); %s

    %Evaporation Equilibration Time (=tau_mix), assume X = 1
    tau_evap(:,:,isize) = tau_mix; %s
    
    %Saturation Concentration associated with tau_evap
    D_air = Dn_air.*(Tmod./T0).^mu1;
    Kn = calc_Knudsen(MW, D_air, Dp(isize), Tmod); %Calculate Knudsen Number
    % Fuchs and Sutugin transition regime correction
    Beta = calc_Beta(Kn, alpha_m);
    csatT = calc_csat_from_tauevap(tau_evap(:,:,isize), Tmod, Dp(isize),...
        D_air, Beta, Xm, MW, rho1, sigma1);

    %Convert csat T to csat(298K)
    cstar(:,:,isize) = csatT./T0.*Tmod./exp(dHvap./R.*(1./T0 - 1./Tmod));
%     cstar(:,:,isize) = csatT;

end
    
% %% Plot Cstar as a function of T and Particle Size for each Water Composition
% for iH2O = 1:length(H2O_range)
%     x = T_range-273.15;  %deg C
%     xticks = [x(1):10:x(end)];
%     xticklbl = cellstr(num2str(xticks(:)));
%     y = log10(Dp); %m
%     yticks = [log10(Dp_lo):1:log10(Dp_hi)];
%     yticklbl = cellstr(num2str(round(yticks(:)), '10^{%d}'));
%     z = log10(squeeze(cstar(:,iH2O,:))');
%     zVec = reshape(z,1,size(z,1)*size(z,2));
%     zticks = [-20:2:20];
%     lzlog = 1;
%     lxlog = 0;
%     lylog = 0;
%     lsave = 1;
%     plot_contour(x, y, z, zticks, lzlog, ...
%         'Temperature (\circC)','Particle Diameter (m)', ...
%         ['Saturation Concentration (@298 K) for x_{H_2O} = ' num2str(H2O_range(iH2O))],...
%         'tex', [x(1),x(end)], [y(1),y(end)], xticks, yticks, xticklbl, yticklbl, ...
%         '../../Figs/Timescale_Analysis/', ...
%         ['cstar_from_TDp_H2O' num2str(iH2O)], [], ...
%         lxlog, lylog, lsave);
% end

%% Plot Cstar as a function of Water and Particle Size for select T
for itemp = 1:length(T_range)
    x = H2O_range;  %mole frac
    xticks = [x(1):10:x(end)];
    xticklbl = cellstr(num2str(xticks(:)));
    y = log10(Dp); %m
    yticks = [log10(Dp_lo):1:log10(Dp_hi)];
    yticklbl = cellstr(num2str(round(yticks(:)), '10^{%d}'));
    z = log10(squeeze(cstar(itemp,:,:))');
    zVec = reshape(z,1,size(z,1)*size(z,2));
    zticks = [-20:2:20];
    lzlog = 1;
    lxlog = 0;
    lylog = 0;
    lsave = 1;
    plot_contour(x, y, z, zticks, lzlog, ...
        'Water Mole Fraction','Particle Diameter (m)', ...
        ['Saturation Concentration (@298 K) for T = ' num2str(T_range(itemp))],...
        'tex', [x(1),x(end)], [y(1),y(end)], xticks, yticks, xticklbl, yticklbl, ...
        '../../Figs/Timescale_Analysis/', ...
        ['cstar_from_TDp_H2O' num2str(itemp)], [], ...
        lxlog, lylog, lsave);
end



