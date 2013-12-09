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
Vars=whos;
PersistentVars=Vars([Vars.persistent]);
PersistentVarNames={PersistentVars.name};
clear(PersistentVarNames{:});
% clear all

%File structure
exp_name = 'CappaMix_01';
DSC_Data_Folder = 'DSC_DATA/';
plotDir = ['../../Figs/DSC/' exp_name '/'];
addpath ../AIOMFAC;


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

%Calculate the equilibrium distribution of each species before the
%Temperature ramp begins
pv_i = find_equilibrium(T_ref, T_i, DSC, MW, dHvap, pstar, Xm_i, X_i);

% Vapor mass concentrations (kg/m3) at initial Temperature and Atm Pressure
cgas(1:nspec) = pv_i.*MW./R./T_i;

% Bulk mass concentrations (kg)
ctot(1:nspec) = X_i .* DSC.mass; %kg %Organic Liquid

% Total concentrations of the species 
% ctot(1:nspec) = cgas + X_i.*c_aer_dist(end);
corg_liq(1:nspec) = ctot - cgas .* DSC.vap_vol;


% Mean velocity of the gas molecules:
c_ave(1:nspec) = sqrt(8.*R.*T_i./MW./pi);
% Diffusion coefficient (m2/s)
D_air(1:nspec) = Dn_air.*(T_i./T_ref).^mu1;
% Mean free path of the gas molecules:
lambda(1:nspec) = 3.*D_air./c_ave;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Masses for All Phases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gas phase concentrations [kg/m3]
Gc(1:nspec) = cgas(1:nspec);
% Organic Luquid total mass [kg]
OLc(1:nspec) = corg_liq(1:nspec);
% Solid Phase total mass [kg]
Sc(1:nspec) = 2.5e-6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculating the time-dependent evaporation                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input = [T_i, Sc, OLc, Gc]'; % T, Temperature, K  
                             % Sc,  Mass of solid phase species (kg)
                             % OLc, Masses of each species bulk (kg)
                             % Gc, Gas phase concentrations (kg/m3)
time = linspace(DSC.time(1), DSC.time(end), 1000);
% time = linspace(DSC.time(1), 2930, 1000);
dt = mean(diff(time));

% Solving the mass fluxes 
options = odeset('RelTol',1E-6,'AbsTol',1E-16);    
[tout, output] = ode45(@fluxes_bulk, time, input,  options, dt,nspec,...
    Cp_liq, Cp_vap,...
    pstar,dHvap,T_ref,MW,sigma1,rhol,Dn_air,mu1,alpha_m, DSC, molec_group_flag,...
    molec_group_stoich);


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
OLc_out = output(1:length(time),nspec+2: 2*nspec+1); %kg
Cp_sys = Cp_sys ./ sum(OLc_out,2);  %dE/dT -> J/(kg degC)

T_out = output(:,1) - 273.15;  %K -> deg C

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Plot output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot Temperature
lylog = 0; %Linear Y Scale
plot_2D({DSC.time, time'}, {DSC.T-273.15, T_out},{1,0},'Time(s)','Temperature(deg C)',...
    'Temperature Profile Comparison',[1,time(end)],[-50,150],...
    plotDir,'DSC_Temp',{'Measured','Modeled'},lylog,1)

%Plot System Heat Capacity
plot_2D({DSC.time, time'}, {DSC.Cp./1000, Cp_sys./1000},{1,0},'Time(s)','Heat Capacity (J/(g degC)',...
    'Heat Capacity Comparison',[1,time(end)],[0,6],...
    plotDir,'DSC_Cp',{'Measured','Modeled'},lylog,1)

%Plot Chemical Composition of All Phases
Sc_out = output(1:length(time),2:nspec+1) .* 1e6; %kg -> mg
Y_limit = ceil(meas.mass/10)*10;
ylog = 0;  %Linear Y-axis
plot_area(time', Sc_out,'Time(s)','Solid Mass (mg)',...
    'Solid Composition Evolution',[1,time(end)],[0,Y_limit],ylog,...
    plotDir,'DSC_Solid_Chem',{},1)


OLc_out = OLc_out .* 1e6; %kg -> mg
Y_limit = ceil(meas.mass/10)*10;
ylog = 0;  %Linear Y-axis
plot_area(time', OLc_out,'Time(s)','Organic Liquid Mass (mg)',...
    'Liquid Composition Evolution',[1,time(end)],[0,Y_limit],ylog,...
    plotDir,'DSC_Liquid_Chem',{},1)

Gc_out = output(1:length(time),2*nspec+2: 3*nspec+1) .* 1e6 .* DSC.vap_vol; %kg/m3 -> mg
Y_limit = max(sum(Gc_out,2))*2;
ylog = 1; %Log Y-axis
plot_area(time', Gc_out,'Time(s)','Mass (mg)',...
    'Gas Chemical Composition Evolution',[1,time(end)],[0,Y_limit],ylog,...
    plotDir,'DSC_Gas_Chem',{},1)

%Recalculate extraneous parameters from ODE function
for ii = 1:length(time)
    [~, p_h2o(ii,:)] = fluxes_bulk(time(ii),output(ii,:)',dt,nspec,Cp_liq,Cp_vap,...
    pstar,dHvap,T_ref,MW,sigma1,rhol,Dn_air,mu1,alpha_m,DSC, molec_group_flag,...
    molec_group_stoich);
end

press = p_h2o(:,1);
h2o_partp = p_h2o(:,2);
h2o_pvap = p_h2o(:,3);
h2o_eqvap = p_h2o(:,4);
h2o_eqliq = p_h2o(:,5);
h2o_flux = p_h2o(:,6);
liq_vol = p_h2o(:,7);
vap_vol = p_h2o(:,8);

%Plot Pressure
lylog = 1; %Log Y-Axis
plot_2D({time', time',time'}, {press, h2o_partp, h2o_pvap},...
    {0,0,0},'Time(s)','Pressure (Pa)',...
    'Pressure Comparison',[1,time(end)],[1e-1,1e6],...
    plotDir,'DSC_Temp',{'Total','P_{H2O}','P_{vap,H2O}'},lylog,1)

%Plot Water Flux
lylog = 0; %Log Y-Axis
plot_2D({time'}, {h2o_flux.*1e6},...
    {0},'Time(s)','Water Flux (mg s^{-1})',...
    'Water Flux',[1,time(end)],[-0.001,0.001],...
    plotDir,'DSC_Temp',{},lylog,1)

%Plot Equlibrium Concentrations
lylog = 0;
plot_2D({time',time',time',time'}, {Gc_out(:,nspec),OLc_out(:,nspec),...
    h2o_eqvap.*1e6,h2o_eqliq.*1e6},{0,0,0,0},'Time(s)','Water Conc (mg)',...
    'Dynamic vs. Equilibrium',[1,time(end)],[0,20],...
    plotDir,'DSC_Temp',{'Gas_{Dyn}','OrgLiq_{Dyn}','Gas_{Eq}','OrgLiq_{Eq}'},lylog,1)

%Print Total Species Evaporation/
for ii = 1:nspec
fprintf(1,'Flow of species %i to liquid phase(mg):  %10.5e\n',ii,OLc_out(end,ii)-OLc_out(1,ii));
end



