function pv_i = find_equilibrium(T_ref, T_i, DSC, MW, dHvap, pstar, Xm_i, X_i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file finds an equilibrium state given current conditions           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% T_ref/T_i: Reference/Initial T [K]
% DSC: structure of calorimeter properties
% MW(1:nspec): molecular weight [kg/mol]
% dHvap(1:nspec): enthalpy of varporization [J/mol]
% pstar(1:nspec): vapor pressure at T_ref
% Xm_i(1:nspec): Mole fraction of each species assuming STP

R = 8.314; %Gas constant [J mol-1 K-1]
nspec = length(MW);

%Don't iterate for the equilibirum initially. Just calculate the vapor
%pressure and subtract from the 

% Equilibrium pressures at the DSC temperature
psat = pstar.*exp(dHvap.*(1./T_ref - 1./T_i)./R);
csat(1:nspec) = psat./R./T_i; %Satn Concentration [mol/m3]

%Total Concentration in the cell [mol species / m3 DSC]
ctot_i(1:nspec) = X_i .* DSC.mass ./ MW ./ DSC.tot_vol;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Vapor-Liquid Equilibrium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use Bisection to find sum(ctot_i) and cliq_i/ctot_i
x1 = 0;
x2 = sum(ctot_i); %Upper Bound on Coa
i = 0;
while (x2-x1) > 0.00001 && i < 100
    xmid = 0.5*(x2+x1);
    psi2 = 1./(1+csat./xmid);
    val = sum(ctot_i .* psi2) - xmid;
    if val < 0
        x2 = xmid;
    else
        x1 = xmid;
    end
    i = i + 1;
end
cliq = ctot_i.*psi2;      %Total moles of liquid [mol/m3]
cgas = ctot_i.*(1-psi2);  %Total moles of vapor [mol/m3]
Xm_i = cliq ./ sum(cliq);  %Mole fraction of each species in the condensed phase


%Convert Back to Vapor Pressures
peq = Xm_i(1:nspec).*psat(1:nspec);

%Equilibrium Gas-Phase Moles of Stuff
moles_vap = DSC.moles_air + sum(cgas) .* DSC.tot_vol;
p_i = moles_vap ./ DSC.vap_vol .* R .* T_i;

pv_i = cgas.*DSC.tot_vol./moles_vap .* p_i;  %Initial Partial pressure at equilibrium
                             %Assume atmospheric pressure first (FIX
                             %THIS!)
                             
                             
end