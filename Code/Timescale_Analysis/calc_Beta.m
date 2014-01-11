function Beta = calc_Beta(Kn,alpha_m)
%Calculate Fuchs-Sutugin Correction factor

Beta = (1.0 + Kn)./...
    (1.0 + (4./alpha_m./3 + 0.377).*Kn + 4.*Kn.^2./alpha_m./3);

end