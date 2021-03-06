function flx0 = fluxes_size_dist(t,input0,dt,nbins,nspec,l_heat,t_res_heat,...
    T_f,T_i,ntot0,n_dist0,pstar,dHvap,T_ref,MW,sigma1,rho,Dn,mu1,press,alpha_m,alpha_t)
    
R = 8.314472;

for i = 1:nspec
    input(1:nbins+3,i) = input0(1+(i-1).*(nbins+3):i.*(nbins+3));
end

% % Accounting for the temperature profile
 xcoord = l_heat.*t./t_res_heat;
% if t_res_heat <17.0 & t_res_heat >13.0 & xcoord < 0.29
% % Temperature profile for 16 s residence time in the beginning of the
% % reactor
%  para = [2.4509 -8.0918 5.0610 0.1405];
%  T_nd = para(1).*xcoord.^3 + para(2).*xcoord.^2 + para(3).*xcoord + para(4);
%  T_TD = T_nd.*(T_f - T_i) + T_i;
% elseif t_res_heat <17.0 & t_res_heat >13.0 & xcoord > 0.46
%  para = [-78.8203 88.5349 -32.7265 4.9935];
%  T_nd = para(1).*xcoord.^3 + para(2).*xcoord.^2 + para(3).*xcoord + para(4);
%  T_TD = T_nd.*(T_f - T_i) + T_i;
% elseif t_res_heat <140 &t_res_heat >100 & xcoord < 0.29
%  para = [-4.0977	-4.0217	4.5377	0.1132];
%  T_nd = para(1).*xcoord.^3 + para(2).*xcoord.^2 + para(3).*xcoord + para(4);
%  T_TD = T_nd.*(T_f - T_i) + T_i;
% elseif t_res_heat <140 &t_res_heat >100 & xcoord > 0.45
%  para = [-72.7449	77.7706	-27.2038	4.1143];
%  T_nd = para(1).*xcoord.^3 + para(2).*xcoord.^2 + para(3).*xcoord + para(4);
%  T_TD = T_nd.*(T_f - T_i) + T_i;
%else    
% Flat temperature profile
T_TD = T_f;
%end

% Accounting for the temperature profile
xcoord0 = l_heat.*(t-dt)./t_res_heat;
% if t_res_heat <17.0 & t_res_heat >13.0 & xcoord0 < 0.29
% % Temperature profile for 16 s residence time in the beginning of the
% % reactor
%  para0 = [2.4509 -8.0918 5.0610 0.1405];
%  T_nd0 = para0(1).*xcoord0.^3 + para0(2).*xcoord0.^2 + para0(3).*xcoord0 + para0(4);
%  T_TD0 = T_nd0.*(T_f - T_i) + T_i;
% elseif t_res_heat <17.0 & t_res_heat >13.0 & xcoord0 > 0.46
%  para0 = [-78.8203 88.5349 -32.7265 4.9935];
%  T_nd0 = para0(1).*xcoord0.^3 + para0(2).*xcoord0.^2 + para0(3).*xcoord0 + para0(4);
%  T_TD0 = T_nd0.*(T_f - T_i) + T_i;
% elseif t_res_heat <140 &t_res_heat >100 & xcoord0 < 0.29
%  para0 = [-4.0977	-4.0217	4.5377	0.1132];
%  T_nd0 = para0(1).*xcoord0.^3 + para0(2).*xcoord0.^2 + para0(3).*xcoord0 + para0(4);
%  T_TD0 = T_nd0.*(T_f - T_i) + T_i;
% elseif t_res_heat <140 &t_res_heat >100 & xcoord0 > 0.45
%  para0 = [-72.7449	77.7706	-27.2038	4.1143];
%  T_nd0 = para0(1).*xcoord0.^3 + para0(2).*xcoord0.^2 + para0(3).*xcoord0 + para0(4);
%  T_TD0 = T_nd0.*(T_f - T_i) + T_i;
% else    
% Flat temperature profile
T_TD0 = T_f;
%end

% Mass of compounds in each particle (kg)
Pc(1:nbins+1,1:nspec) = input(1:nbins+1,1:nspec);
% Gas phase mass concentrations at the TD temperature
Gc(1,1:nspec) = input(nbins+2,1:nspec).*T_TD0./T_TD;
Gc(2,1:nspec) = input(nbins+3,1:nspec).*T_TD0./T_TD;
% Total aerosol number concentration (1/m3)
ntot = ntot0.*T_i./T_TD;
% Aerosol numbers in each of the bins (1/m3)
n_dist = n_dist0.*T_i./T_TD;

for i = 1:nspec
    for j = 1:nbins+1
    if Pc(j,i) <= 0.0 || ~isreal(Pc(j,i))
    Pc(j,i) = 0;
    end
    if Gc(i) <= 0.0 || ~isreal(Gc(i))
    Gc(i) = 0;
    end
    end
end

% Equilibrium pressures at the TD temperature
psat(1:nspec) = pstar.*exp(dHvap.*(1./T_ref - 1./T_TD)./R);
csat(1:nspec) = MW.*psat./R./T_TD;

% Total particle mass at each size bin
mp(1:nbins+1,1) = sum(Pc,2);

% Diffusion coefficients of the species
D(1:nspec) = Dn(1:nspec).*(T_TD./T_ref).^mu1(1:nspec);

% Calculating the composition
for j = 1:nbins+1
    
    if mp(j) > 1.0e-22;
        % Mass fractions of each species in the particles
        X(j,1:nspec) = Pc(j,1:nspec)./mp(j);
        % Mole fractions
        n_apu(1:nspec) = X(j,1:nspec)./MW(1:nspec);
        n_tot_apu = sum(X(j,1:nspec)./MW(1:nspec));
        Xm(j,1:nspec) = n_apu(1:nspec)./n_tot_apu;
        % Particle density
        rhol(j) = sum(X(j,1:nspec).*rho(1:nspec));
        % Particle surface tension
        sigmal(j) = sum(Xm(j,1:nspec).*sigma1(1:nspec));
        % Particle volume and size
        vp(j) = mp(j)./rhol(j);
        rp(j) = (3.*vp(j)./4./pi).^(1./3);
        % Kelvin effect
        Ke(j,1:nspec) = exp(2.0.*MW(1:nspec).*sigmal(j)./R./T_TD./rhol(j)./rp(j));
        % Equilibrium pressures, mole based
        %peq(j,1:nspec) = Xm(j,1:nspec).*psat(1:nspec).*Ke(j,1:nspec);
        % Mass based! Remember to fix back!!
        peq(j,1:nspec) = X(j,1:nspec).*psat(1:nspec).*Ke(j,1:nspec);
        
        
        % Transitional correction
        % Mean velocity of the gas molecules:
        c_ave(1:nspec) = sqrt(8.*R.*T_TD./MW(1:nspec)./pi);
        % Mean free path of the gas molecules:
        lambda(1:nspec) = 3.*D(1:nspec)./c_ave(1:nspec);
        % Knudsen number:
        Kn(j,1:nspec) = lambda(1:nspec)./rp(j);
        % Fuchs and Sutugin transition regime correction
        beeta(j,1:nspec) = (1.0 + Kn(j,1:nspec))./(1.0 + (4./alpha_m(1:nspec)./3 + 0.377).*Kn(j,1:nspec) + 4.*Kn(j,1:nspec).^2./alpha_m(1:nspec)./3);

    else
        mp(j) = 0.0;
        vp(j) = 0.0;
        rp(j) = 0.0;
        peq(j,1:nspec) = zeros(1,nspec);
        beeta(j,1:nspec) = zeros(1,nspec);
    end
end

% Pressure at the particle surface
pv_a(1:nbins+1,1:nspec) = peq;

% Partial pressures far away from the particles
pv_i(1,1:nspec) = Gc(1,1:nspec).*R.*T_TD./MW(1:nspec);


for j = 1:nbins+1
    for i = 1:nspec
        if mp(j) > 0.0 & (1.0 - pv_a(j,i)./press)/(1.0 - pv_i(1,i)./press) > 0.0
            % Mass flux of each compound to the particles
            flx(j,i) = 4.*pi.*rp(j).*press.*D(i).*beeta(j,i).*MW(i).*log((1.0 - pv_a(j,i)./press)./(1.0 - pv_i(1,i)./press))./R./T_TD;
        elseif mp > 0.0 & (1.0 - pv_a(j,i)./press)/(1.0 - pv_i(1,i)./press) <= 0.0
            flx(j,i) = 4.*pi.*rp(j)*D(i).*beeta(j,i).*MW(i).*(pv_i(1,i)-pv_a(j,i))./R./T_TD;
        else
            flx(j,i) = 0.0;
        end
    end
    % flx_apu_dist(j,1:nspec) = n_dist(j,1).*flx(j,1:nspec);
    flx_apu_dist(j,1:nspec) = n_dist(j,1).*flx(j,1:nspec);
end

% Mass flux of each compound to the gas phase
flx(nbins+2,1:nspec) = -sum(flx_apu_dist(1:end-1,1:nspec),1);
flx(nbins+3,1:nspec) = -sum(ntot.*flx(nbins+1,1:nspec));

% flx(1:nbins+1,1:9)
% flx(nbins+2:nbins+3,1:9)


% Changing to column vector
for i = 1:nspec
    flx0(1+(i-1)*(nbins+3):i*(nbins+3)) = flx(1:nbins+3,i);
end

flx0 = flx0';
for ll = 1:length(flx0)
    if ~isreal(flx0(ll)) | isnan(flx0(ll))
        flx0(ll)=0.0;
    end
end


end