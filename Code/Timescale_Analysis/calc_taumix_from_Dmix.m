function tau_mix = calc_taumix_from_Dmix(Dp, Dmix)
%Calculate Mixing Timescale given Particle diameter and Diffusivity
% tau = r^2/ (pi^2 * Dmix);

pi = 3.1415;
tau_mix = (Dp./2).^2 ./ (pi.^2) ./ Dmix;

end



%X = tau_evap / tau_mix

%X = rho1 .* pi .* Dmix ./
%    (Beta .* 4.* D_air .* csat .* Xm .* Ke)
