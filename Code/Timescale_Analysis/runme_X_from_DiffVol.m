%Script for creating plots of X from values of Diffusivity and Volatility
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
Dmix_lo = 1.0e-21;  %m2 s-1
Dmix_hi = 1.0e-9;   %m2 s-1
% Dmix = logspace(log10(Dmix_lo),log10(Dmix_hi),50);  %m2 s-1
Dmix = 10.^[-20 -18 -15 -11 -10];  %m2 s-1

%Volatility Range
Csat_lo = 1.0e-8;  %ug m-3
Csat_hi = 1.0e7;   %ug m-3
Csat = logspace(log10(Csat_lo),log10(Csat_hi),60);  %ug m-3

%Size Range
Dp_lo = 1.0e-9; %m
Dp_hi = 10.0e-6; %m
Dp = logspace(log10(Dp_lo),log10(Dp_hi),50);  %m
% Dp = 10.^[-9 -8 -7 -6 -5]; %m

%Loop through several sizes
for isize = 1:length(Dp)


    %Mixing Equilibration Time (fxn of Dmix and Size)
    tau_mix = calc_taumix_from_Dmix(Dp(isize),Dmix); %s

    %Evaporation Equilibration Time (fxn of Vol and Size)
    T = T0;
    D_air = Dn_air.*(T/T0).^mu1;
    Kn = calc_Knudsen(MW, D_air, Dp(isize), T); %Calculate Knudsen Number
    % Fuchs and Sutugin transition regime correction
    Beta = calc_Beta(Kn, alpha_m);
    tau_evap = calc_tauevap_from_csat(Csat, T, Dp(isize), D_air, Beta, Xm, MW,...
        rho1, sigma1); %s

    %Calculate X = tau_evap / tau_mix
    Xi(:,:,isize) = tau_evap' * (1./tau_mix);

end    


% %Plot X as a function of Diffusivity and Volatility
% for isize = 1:length(Dp)
    % x = Dmix .* 100^2; %cm2 s-1
    % xticks = 10.^[log10(Dmix_lo.*100^2):1:log10(Dmix_hi.*100^2)];
    % y = Csat;  %ug m-3
    % yticks = 10.^[-5:2:8];
    % z = log10(Xi(:,:,isize));
    % zVec = reshape(z,1,size(z,1)*size(z,2));
    % zticks = [-5:3:20];
    % zlbl = zticks;
    % lxlog = 1;
    % lylog = 1;
    % lsave = 1;
    % plot_contour(x, y, z, zticks, zlbl, 'Diffusivity (cm^2 s^{-1})', ...
    %     'Saturation Concentration (\mu g m^{-3})', ...
    %     ['X (= \tau_{evap}/\tau_{mix}) for particles of size ' num2str(Dp(isize)*1.0e6) ' \mum'],...
    %     'tex', [x(1),x(end)], [y(1),y(end)], xticks, yticks, 'Figs/', ...
    %     ['X_from_DiffVol_Dp' num2str(isize)], [], ...
    %     lxlog, lylog, lsave);
% end

%Plot X as a function of Size and Volatility
for imix = 1:length(Dmix)
    x = Csat;  %ug m-3
    xticks = 10.^[log10(Csat_lo):2:log10(Csat_hi)];
    y = Dp; %m
    yticks = 10.^[log10(Dp_lo):1:log10(Dp_hi)];
    z = log10(squeeze(Xi(:,imix,:))');
    zVec = reshape(z,1,size(z,1)*size(z,2));
    zticks = [-10:2:20];
    zlbl = zticks;
    lxlog = 1;
    lylog = 1;
    lsave = 1;
    plot_contour(x, y, z, zticks, zlbl, ...
        'Saturation Concentration (\mu g m^{-3})','Size (m)', ...
        ['X (= \tau_{evap}/\tau_{mix}) for Diffusivity ' num2str(Dmix(imix)*100^2) ' cm^2 s^{-1}'],...
        'tex', [x(1),x(end)], [y(1),y(end)], xticks, yticks, 'Figs/', ...
        ['X_from_DpVol_Diff' num2str(imix)], [], ...
        lxlog, lylog, lsave);
end






