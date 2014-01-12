function tau_evap = calc_tauevap_from_csat(csat, T, Dp, D_air, Beta, Xm, ...
    MW, rho1, sigma1)
%Calculates equilibrium saturation concentration given an evaporation
%timescale
%
% csat - Saturation concentration (ug m-3)
% T - temperature (K)
% Dp - particle diameter (m)
% D_air - diffusivity in air (m^2/s)
% Beta - Fuchs and Sutugin correction factor
% Xm - mole-fraction of species in the particle phase
% MW - molecular weight (kg mol-1)
% rho - condensed-phase density (kg cm-3)

R = 8.314;  %J mol-1 K-1
rp = Dp ./ 2;  %particle radius

%Saturation Pressure (Pa) at T
psat = csat ./MW.*R.*T./1.0e9;  %Pa

%Equilibrium Pressure (Pa) for Size rp
Ke = exp(2.0.*MW.*sigma1./R./T./rho1./rp);  % Kelvin effect
peq = psat .* Xm .* Ke;

%Calculate Equilibrium Vapor Pressure
% tau_evap = R.*T.*rho1 .* rp.^2 ./ ...
%     (Beta .* 4.*pi.* MW .* D_air .* peq); %s
tau_evap = R.*T ./ ...
    (Beta .* 4.*pi.* MW .* D_air .* rp .* peq); %s


end