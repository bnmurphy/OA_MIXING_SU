%This script imports processed MFRs from the trajectory model and 
%observations and makes plots of them

age = {'FULL_HYSP'};
% tags = {'2bin_nohet_72';'FUNC100_soft_nohet_72'};
% ctags = {'2bin-nohet-72'; 'FUNC100-soft-nohet-72'};
tags = {'2bin_ksivoc_1stgen'};
ctags = {'2bin-ksivoc-1stgen'};

dhvaps = [50 75 100];
accom = [0.01 0.1 1];
days = {'avg'};

ntags = size(tags,2);
ndh = length(dhvaps);
naccom = length(accom);
ndays = size(days,2);
nhr = 24;

%Get All Data
mfr = zeros(ntags, ndh, naccom, nhr, ndays);
get_mfr
set_meas


%Make Plots of the Data
make_plots_td


