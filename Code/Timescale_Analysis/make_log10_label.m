function log10_label = make_log10_label( exps )
%Concatenate labels of powers of 10 to use for plot axes

ntick = length(exps);

log10_label = cell(1,ntick);
for itick = 1:ntick
    log10_label(itick) = {['10^{' num2str(exps(itick)) '}']};
end

end