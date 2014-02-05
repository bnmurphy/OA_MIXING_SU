function [pv_i, mliq_i, msld_i, Xm_i, X_i] = find_equilibrium(Tmelt, T_ref, T_i, DSC, MW, dHfus, dHvap,...
    pstar, Xm_i, X_i, molec_group_flag, molec_group_stoich)
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

%Get Activity Coefficients
% ln_gamma = AIOMFAC_gamma_SR_v1(nspec, Xm_i, molec_group_flag, molec_group_stoich, T_i);
% gamma = exp(nonzeros(ln_gamma))';  %Convert activity coeff sparse matrix to vector  
gamma = ones(1,nspec);

% Equilibrium pressures at the DSC temperature
psat = gamma .* pstar.*exp(dHvap.*(1./T_ref - 1./T_i)./R);
csat(1:nspec) =  psat./R./T_i; %Satn Concentration [mol/m3]

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
cliq_i = ctot_i.*psi2;      %Total moles of liquid [mol/m3]
cgas_i = ctot_i.*(1-psi2);  %Total moles of vapor [mol/m3]
Xm_i = cliq_i ./ sum(cliq_i);  %Mole fraction of each species in the liquid phase

%Convert Back to Vapor Pressures
peq = Xm_i(1:nspec).*psat(1:nspec);

%Equilibrium Gas-Phase Moles of Stuff
moles_vap = DSC.moles_air + sum(cgas_i) .* DSC.tot_vol;
p_i = moles_vap ./ DSC.vap_vol .* R .* T_i; %Total Gas Pressure Initially

pv_i = cgas_i.*DSC.tot_vol./moles_vap .* p_i;  %Initial Partial pressure at equilibrium
                             

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Solid-Liquid Equilibrium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calcuate solubility (xsol) as a fxn of T_i (moles solute/mole solvent)
xsol_i(1:nspec) = 1./gamma(1:nspec) .* exp(dHfus./R .* (1./Tmelt - 1./T_i));
cliq_i0 = cliq_i;
%Use Bisection to find:
%   cliq_tot - total liquid concentration
%   ctot     - total liquid+solid conc of all species
%   ctot_i   - total liquid+solid of each species
%   csld_i   - solid-phase concentration of each species
guess1 = 0;  %Lower bound on total liquid guess
guess2 = sum(ctot_i); %Upper Bound on liquid guess
i = 0;
while (guess2-guess1) > 0.00001 && i < 100
    guess_mid = 0.5*(guess2+guess1);
    cliq_i = min(cliq_i0, xsol_i .* guess_mid);
    val = sum(cliq_i) - guess_mid;
    if val < 0
        guess2 = guess_mid;
    else
        guess1 = guess_mid;
    end
    i = i + 1;
end
csld_i = ctot_i - cliq_i; %Total moles of solid [mol/m3]
%cliq_i = ctot_i - csld_i; %Total moles of liquid [mol/m3]

Xm_i = cliq_i ./ sum(cliq_i); %Mole fraction of each species in the liquid phase

%Calculate Liquid Phase Mass Fractions
X_i = Xm_i.*MW ./ sum( Xm_i.*MW );

mliq_i = cliq_i .* MW .* DSC.liq_vol; %kg of liquid
msld_i = csld_i .* MW .* DSC.liq_vol; %kg of solid

                            
end