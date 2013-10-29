%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This file contains the properties of the initial size distribution                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Size bins
% Number of bins
nbins = 5;

% Size range (m)
% dprange = [50e-9];

% Bin Diameters (m)
% dp_int = logspace(log10(dprange(1)),log10(dprange(end)),nbins)
dp_int = [0.0586, 0.1172, 0.2344, 0.4688, 0.9375] .* 1e-6;
rp_int = dp_int./2;

% Masses in each bin (kg/m3)
% c_aer_dist_int = c_aer_tot_i;
% c_aer_dist_int = [a b c d e ];

% Number in each bin (1/m3)
% n_dist_int = n_tot_i;
n_dist_int = c_aer_dist_int./(4.*1400.*rp_int.^3.*pi./3);

% Total mass (kg/m3)
c_aer_int_i = nansum(c_aer_dist_int);
n_int_i = nansum(n_dist_int);

% Monodisperse case total aerosol mass concentration (kg/m3)
% c_aer_tot_i = 13.0./1e9;
c_aer_tot_i = nansum(c_aer_dist_int);

% Total aerosol number (1/m3)
%n_tot_i = 4e5.*1e6;

% Peak size (m)
% dp_peak_i = 50e-9;
dp_peak_i = dp_int( find( c_aer_dist_int == max(c_aer_dist_int) ) );
% dp_peak_i = dp_int( find( n_dist_int == max(n_dist_int) ) );

rp_peak_i = dp_peak_i./2;
% Here note that the density is assumed to be 1400 kg/m3
n_tot_i = c_aer_tot_i./(4.*1400.*rp_peak_i.^3.*pi./3);

% Collecting all to one vector 
dp_i = [dp_int dp_peak_i]';
rp_i = dp_i./2;
% c_aer_dist_i = [c_aer_dist_int c_aer_tot_i]';
% n_dist_i = [n_dist_int n_tot_i]';
c_aer_dist_i = [c_aer_dist_int c_aer_int_i]';
n_dist_i = [n_dist_int n_int_i]';







