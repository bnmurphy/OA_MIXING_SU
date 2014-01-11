function Kn = calc_Knudsen(MW, D_air, rp, T)
%Calculate Knudsen Number given molecular weight and particle radius
% MW - molecular weight kg mol-1
% D_air - Diffusivity in air m2/s
% rp - particle radius m
% T - temperature K

pi = 3.1415;
R = 8.314;

c_ave = sqrt(8.*R.*T./MW./pi); % Mean velocity of the gas molecules
lambda = 3.*D_air./c_ave;      % Mean free path of the gas molecules
Kn = lambda./rp;          % Knudsen number

end