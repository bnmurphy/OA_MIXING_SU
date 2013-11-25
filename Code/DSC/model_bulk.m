%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main program to model the dynamic evaporation in the        %
%   Differential Scanning Calorimeter. An energy profile and mixture      %
%   properties are input to the model and temperature and concentration   %
%   profiles are output
%
% Although the DSC experiment focuses on bulk fluid, try to keep particle
% functionality throughout the code, so that it can be extended to look at
% particle issues afterward
%
% Experiment Names:
%     CappaMix_01
%     CappaMix_02
%     CappaMix_03_H2O_025
%     CappaMix_04_H2O_025
%     CappaMix_06_H2O_050
%     CappaMix_07_H2O_075
%     CappaMix_08_H2O_075
%     CappaMix_09_H2O_090
%     CappaMix_10_H2O_090
%     CappaMix_11_H2O_095
%     CappaMix_12_H2O_095
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear variables
% clear all

%File structure
exp_name = 'CappaMix_12_H2O_095';
DSC_Data_Folder = 'DSC_DATA/';
plotDir = ['../../Figs/DSC/' exp_name '/'];

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
    % Mole fraction based equilibrium pressures
    %peq_i(1:nbins+1,i) = Xm_i(i).*psat_i(i).*Ke_i(1:nbins+1,i);
    
    % Mass fraction based equilibrium pressures
    peq_i(i) = X_i(i).*psat_i(i);
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

% Bulk mass concentration (kg)
Bc = mp_i;


% Total concentrations of the species 
% ctot(1:nspec) = cgas + X_i.*c_aer_dist(end);
ctot(1:nspec) = cgas + X_i .* Bc;
    
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
time = linspace(DSC.time(1), DSC.time(end), 1000);
dt = mean(diff(time));

% Solving the mass fluxes 
options = odeset('RelTol',1E-6,'AbsTol',1E-10);    
[tout, output] = ode45(@fluxes_bulk, time, input,  options, dt,nspec,...
    Cp_liq, Cp_vap,...
    pstar,dHvap,T_ref,MW,sigma1,rho,Dn,mu1,p,alpha_m, DSC);


%Compute experimental values
%First get Temperature from model
T_out = output(:,1); %K
%Find change in Temperature with time
dT = diffxy(time, T_out);  %dT/dt -> degC/s

%Retrieve Heat Flow from Experiment
Q = zeros(length(time),1);
for itime = 1:length(time)
   Q(itime) = mean(DSC.Q(find(DSC.time < time(itime),1,'last' ): ...
                          find(DSC.time > time(itime),1,'first')  ) );  %J/s
   if isnan(Q(itime))
       if time(itime) <= DSC.time(1), Q(itime) = DSC.Q(1); end
       if time(itime) >= DSC.time(end), Q(itime) = DSC.Q(end); end
   end
end

%Divide Heat Flow by Temperature Temporal Evolution
%to get heat capacity
Cp_sys = Q ./ dT; %dE/dT -> J/degC

%Normalize with experimental mass to get specific heat capacity
Cp_sys = Cp_sys ./ DSC.mass;  %dE/dT -> J/(kg degC)


T_out = output(:,1) - 273.15;  %K -> deg C

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Plot output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot Temperature
plot_2D({DSC.time, time'}, {DSC.T-273.15, T_out},{1,0},'Time(s)','Temperature(deg C)',...
    'Temperature Profile Comparison',[1,time(end)],[-50,150],...
    plotDir,'DSC_Temp',{'Measured','Modeled'},1)

%Plot System Heat Capacity
plot_2D({DSC.time, time'}, {DSC.Cp./1000, Cp_sys./1000},{1,0},'Time(s)','Heat Capacity (J/(g degC)',...
    'Heat Capacity Comparison',[1,time(end)],[0,6],...
    plotDir,'DSC_Cp',{'Measured','Modeled'},1)

%Plot Chemical Composition of Both Phases
Bc_out = output(1:length(time),2:nspec+1) .* 1e6; %kg -> mg
Y_limit = ceil(meas.mass/10)*10;
ylog = 0;  %Linear Y-axis
plot_area(time', Bc_out,'Time(s)','Mass (mg)',...
    'Bulk Chemical Composition Evolution',[1,time(end)],[0,Y_limit],ylog,...
    plotDir,'DSC_Cond_Chem',{},1)

Gc_out = output(1:length(time),nspec+2: 2*nspec+1) .* 1e6 .* DSC.vap_vol; %kg/m3 -> ug/m3
Y_limit = max(sum(Gc_out,2))*2;
ylog = 1; %Log Y-axis
plot_area(time', Gc_out,'Time(s)','Mass (mg)',...
    'Gas Chemical Composition Evolution',[1,time(end)],[0,Y_limit],ylog,...
    plotDir,'DSC_Gas_Chem',{},1)

% %Add Temperature Axis
% hold on 
% axT = axes;
% h = plot(time',T_out,'linewidth',2);
% set(axT,'ylim',[-50,150],'xlim',[1,time(end)]);

b = ou2;




