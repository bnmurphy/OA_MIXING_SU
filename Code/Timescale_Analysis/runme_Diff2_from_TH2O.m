%Script to plot the Diffusivity of condensed phase as a function of T and
%H2O content

clear variables

%% Setup Physical Parameters
pi = 3.1415;
R = 8.314; %J mol-1 K-1
k = 1.38e-23; %m2 kg s-2 K-1
Dn_air = 5.0e-6;  %m2 s-1
mu1 = 1.75;  %Parameter for D_air calculation
rh = 2.12e-10; %Hydrodynamic Radius (m)

%% Load Viscosity Measurements
load('expt_visc.mat')


%% Set Independent Variable Ranges
%Temperature Range
T_lo = 273.15; %K
T_hi = 90+273.15; %K
T_range = linspace(T_lo,T_hi,900);

%Water Composition Range
H2O_lo = 0.0;
H2O_hi = 0.95;
H2O_range = linspace(0.0, 0.95, 1000);

%% Find Diffusivity - interpolate H2O and T
H2Omod = repmat(H2O_range,length(T_range),1);
Tmod = repmat(T_range',1,length(H2O_range));
visc.adjust = interp2(expt.T'+273.15, expt.water, expt.visc.adjust', Tmod, H2Omod,'linear');

%Smooth the data
for icol = 1:length(H2O_range)
    visc.smooth(:,icol) = 10 .^ smooth(log10(visc.adjust(:,icol)),'moving');
end
for irow = 1:length(T_range)
    visc.smooth(irow,:) = 10 .^ smooth(log10(visc.smooth(irow,:)),'moving');
end

%% Plot Viscosity
x = T_range-273.15;  %deg C
xticks = [x(1):10:x(end)];
y = H2O_range; %m
yticks = [0.0:0.2:1.0];
z = log10(visc.smooth');
zVec = reshape(z,1,size(z,1)*size(z,2));
zticks = [-1:1:7];
zlbl = zticks;
lxlog = 0;
lylog = 0;
lsave = 1;
plot_contour(x, y, z, zticks, zlbl, ...
    'Temperature (\circC)','Water Mole Fraction', ...
    'Viscosity of Condensed Phase (Pa s)',...
    'tex', [x(1),x(end)], [y(1),y(end)], xticks, yticks, 'Figs/', ...
    'Visc_from_TH2O', [], ...
    lxlog, lylog, lsave);


%Calculate Diffusivity
% lH2O = length(expt.water); 
% lT = length(expt.T);
% expt.Tmat = repmat(expt.T',1,lH2O) + 273.15;  %K
% expt.Dmix = k .* expt.Tmat ./ 6 ./ pi ./ visc.smooth ./ rh; %m2 s-1
Dmix = k .* Tmod ./ 6 ./ pi ./ visc.smooth ./ rh; %m2 s-1

%% Plot Diffusivity
x = T_range-273.15;  %deg C
xticks = [x(1):10:x(end)];
y = H2O_range; %m
yticks = [0.0:0.2:1.0];
z = log10(Dmix'.*100^2);
zVec = reshape(z,1,size(z,1)*size(z,2));
zticks = [-22:1:-6];
zlbl = zticks;
lxlog = 0;
lylog = 0;
lsave = 1;
plot_contour(x, y, z, zticks, zlbl, ...
    'Temperature (\circC)','Water Mole Fraction', ...
    'Diffusivity of Condensed Phase (cm^{2} s^{-1})',...
    'tex', [x(1),x(end)], [y(1),y(end)], xticks, yticks, 'Figs/', ...
    'Dmix_from_TH2O', [], ...
    lxlog, lylog, lsave);