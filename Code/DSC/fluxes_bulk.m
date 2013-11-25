function flx0 = fluxes_bulk(t,input,dt,nspec,Cp_liq,Cp_vap,...
    pstar,dHvap,T_ref,MW,sigma1,rho,Dn,mu1,press,alpha_m,DSC)

if mod(t,100) == 0, t, end

R = 8.314472;
%Load Local Variables
T = input(1);  %Kelvin
Bc = input(2:nspec+1)';  %Bulk concentration kg
Gc = input(nspec+2:2*nspec+1)';   %Gas concentration kg m-3
Qin = DSC.Q(find(DSC.time >= t,1 ) );  %J/s



Bc( Bc <= 0.0 | ~isreal(Bc) ) = 0.0;
Gc( Gc <= 0.0 | ~isreal(Gc) ) = 0.0;

% Equilibrium pressures at the DSC temperature (should be Pa)
psat(1:nspec) = pstar.*exp(dHvap.*(1./T_ref - 1./T)./R);
csat(1:nspec) = MW.*psat./R./T;

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
    
    % Bulk density (Mass-weighted)
    rhol = sum(X(1:nspec).*rho(1:nspec));
    
    % Bulk surface tension (Mole-weighted)
    sigmal = sum(Xm(1:nspec).*sigma1(1:nspec));
    
    % Bulk volume and size
    vp = Bc_tot ./ rhol;
    rp = (3.*vp./4./pi).^(1./3);
    
    % Equilibrium pressures, mole based
    %peq(j,1:nspec) = Xm(j,1:nspec).*psat(1:nspec);
    % Mass based! Remember to fix back!!
    peq(1:nspec) = X(1:nspec).*psat(1:nspec);
    
else
    mp = 0.0;
    vp = 0.0;
    rp = 0.0;
    peq(1:nspec) = zeros(1,nspec);
    beeta(1:nspec) = zeros(1,nspec);
end

% Pressure at the bulk surface (Pa)
pv_a(1:nspec) = peq;

% Partial pressures far away from the bulk (Pa)
pv_i(1:nspec) = Gc(1:nspec).*R.*T./MW(1:nspec);

% Mass flux of each compound to the particles   kg s-1
b_flx = zeros(1,nspec);
for i = 1:nspec
    ln_fact = (1.0 - pv_a(i)./press)/(1.0 - pv_i(i)./press);
    if Bc_tot > 0.0 
        
        if ln_fact > 0.0  %Evaporation
                          %Stefan (Convective) Flow + Diffusive Flow
            
            b_flx(i) = DSC.L.*press.*D(i).*MW(i).* ...
                log((1.0 - pv_a(i)./press)./(1.0 - pv_i(1,i)./press))./R./T;
            
        else        %Condensation
                    %Just Diffusive flow
            
            b_flx(i) = DSC.L.*D(i).*MW(i).* ...
                (pv_i(i)-pv_a(i))./R./T;
        end
        
    else    %Not enough mass; Do Nothing
        b_flx(i) = 0.0;
    end
end

% Mass flux of each compound to the gas phase
g_flx(1:nspec) = -b_flx ./ DSC.tot_vol; %kg m-3 s-1

if T > 373
    t;
end

%Energy Balance (Q is going into the system)
dT_latent = sum( -b_flx(1:nspec) .* dHvap(1:nspec) ./ MW(1:nspec) ); %J s-1
dT_Cpliq = sum( Bc(1:nspec) .* Cp_liq(1:nspec) ./ MW(1:nspec) ); %J K-1
dT_Cpvap = sum( Gc(1:nspec) .* DSC.vap_vol .* Cp_vap(1:nspec)./ MW(1:nspec)); %J K-1
dTdt = ( Qin - dT_latent) / (dT_Cpliq + dT_Cpvap);  %K s-1

if dTdt ~= dTdt
    pause
end

% Changing to column vector
flx0 = [dTdt, b_flx, g_flx]';
flx0( ~isreal(flx0) | isnan(flx0) ) = 0.0;




end