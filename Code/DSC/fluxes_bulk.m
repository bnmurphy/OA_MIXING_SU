function flx0 = fluxes_bulk(t,input,dt,nspec,...
    pstar,dHvap,T_ref,MW,sigma1,rho,Dn,mu1,press,alpha_m,DSC.Q)
    
R = 8.314472;

%Load Local Variables
T = input(1);
Bc = input(2:nspec+1);
Gc = input(nspec+2:2*nspec+1);

%Set power
if t < DSC.time(1)
    Q = 0.0; %No Power before experiment
else
    Q = DSC.Q(find(DSC.time >= t,1 ) ); 
end


Bc( Bc <= 0.0 || ~isreal(Bc) ) = 0.0;
Gc( Gc <= 0.0 || ~isreal(Gc) ) = 0.0;

% Equilibrium pressures at the DSC temperature
psat(1:nspec) = pstar.*exp(dHvap.*(1./T_ref - 1./T_TD)./R);
csat(1:nspec) = MW.*psat./R./T_TD;

% Diffusion coefficients of the species
D(1:nspec) = Dn(1:nspec).*(T/T_ref).^mu1(1:nspec);

% Calculating the composition    
Bc_tot = sum(Bc);
if sum(Bc) > 1.0e-22;
    % Mass fractions of each species in the particles
    X(1:nspec) = Bc(1:nspec) ./ Bc_tot;
    
    % Mole fractions
    n_apu(1:nspec) = X(1:nspec) ./ MW(1:nspec);
    n_tot_apu = sum(X(1:nspec) ./ MW(1:nspec));
    Xm(1:nspec) = n_apu(1:nspec) ./ n_tot_apu;
    
    % Particle density (Mass-weighted)
    rhol = sum(X(1:nspec).*rho(1:nspec));
    
    % Particle surface tension (Mole-weighted)
    sigmal = sum(Xm(1:nspec).*sigma1(1:nspec));
    
    % Particle volume and size
    vp = Bc_tot ./ rhol;
    rp = (3.*vp./4./pi).^(1./3);
    
    % Kelvin effect (nonexistent)
%     Ke(1:nspec) = exp(2.0.*MW(1:nspec).*sigmal./R./T_TD./rhol./rp);
    Ke(1:nspec) = 1.0;
    
    % Equilibrium pressures, mole based
    %peq(j,1:nspec) = Xm(j,1:nspec).*psat(1:nspec).*Ke(j,1:nspec);
    % Mass based! Remember to fix back!!
    peq(1:nspec) = X(1:nspec).*psat(1:nspec).*Ke(1:nspec);
    
    
    % Transitional correction
    % Mean velocity of the gas molecules:
    c_ave(1:nspec) = sqrt(8.*R.*T_TD./MW(1:nspec)./pi);
    % Mean free path of the gas molecules:
    lambda(1:nspec) = 3.*D(1:nspec)./c_ave(1:nspec);
    % Knudsen number:
    Kn(1:nspec) = lambda(1:nspec)./rp;
    % Fuchs and Sutugin transition regime correction
    beeta(1:nspec) = (1.0 + Kn(1:nspec))./(1.0 + (4./alpha_m(1:nspec)./3 + 0.377).*Kn(1:nspec) + 4.*Kn(1:nspec).^2./alpha_m(1:nspec)./3);
    
else
    mp = 0.0;
    vp = 0.0;
    rp = 0.0;
    peq(1:nspec) = zeros(1,nspec);
    beeta(1:nspec) = zeros(1,nspec);
end

% Pressure at the bulk surface
pv_a(1:nspec) = peq;

% Partial pressures far away from the bulk
pv_i(1:nspec) = Gc(1:nspec).*R.*T./MW(1:nspec);


for i = 1:nspec
    if mp > 0.0 && (1.0 - pv_a(i)./press)/(1.0 - pv_i(i)./press) > 0.0
        % Mass flux of each compound to the particles
        b_flx(i) = 4.*pi.*rp.*press.*D(i).*beeta(i).*MW(i).*log((1.0 - pv_a(i)./press)./(1.0 - pv_i(1,i)./press))./R./T;
    elseif mp > 0.0 && (1.0 - pv_a(i)./press)/(1.0 - pv_i(i)./press) <= 0.0
        b_flx(i) = 4.*pi.*rp*D(i).*beeta(i).*MW(i).*(pv_i(i)-pv_a(i))./R./T;
    else
        b_flx(i) = 0.0;
    end
end

% Mass flux of each compound to the gas phase
g_flx(1:nspec) = -b_flx;

%Energy Balance (Q is going out of the system)
dT_latent = sum( b_flx(1:nspec) .* dHvap(1:nspec) );
dT_Cpliq = sum( Bc(1:nspec) .* Cp_liq(1:nspec));
dT_Cpvap = sum( Gc(1:nspec) .* Cp_vap(1:nspec));
dTdt = ( -Q - dT_latent) / (dT_Cpliq + dT_Cpvap);

% Changing to column vector
flx0 = [dTdt, b_flx, g_flx]';
flx0( ~isreal(flx0) || isnan(flx0) ) = 0.0;




end