function csat = calc_csat_from_tauevap(tau_evap, T, rp, D_air, Beta, Xm, ...
    MW, rho1, sigma1)
%Calculates equilibrium saturation concentration given an evaporation
%timescale
%
% tau_evap - evaporation timescale (s)
% T - temperature (K)
% rp - particle radius (m)
% D_air - diffusivity in air (m^2/s)
% Beta - Fuchs and Sutugin correction factor
% Xm - mole-fraction of species in the particle phase
% MW - molecular weight (kg mol-1)
% rho - condensed-phase density (kg cm-3)

R = 8.314;  %J mol-1 K-1

%Calculate Equilibrium Vapor Pressure
peq = R.*T.*rho1 .* rp.^2 ./ ...
    (Beta .* 4.*pi.* MW .* D_air .* tau_evap); %Pa

Ke = exp(2.0.*MW.*sigma1./R./T./rho1./rp);  % Kelvin effect

% Saturation Pressure (Pa)
psat = peq ./ Xm ./ Ke;

%Saturation Concentration (ug m-3) at T0
csat = MW.*psat./R./T .* 1.0e9;

end