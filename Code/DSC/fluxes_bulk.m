function [flx0, p_h2o] = fluxes_bulk(t,input,dt,nspec,Cp_liq,Cp_vap,...
    pstar,dHvap,T_ref,MW,sigma1,rhol,Dn_air,mu1,alpha_m,DSC)

R = 8.314472;
%Load Local Variables
T = input(1);  %Kelvin
OLc = input(2:nspec+1)';  %Organic Liquid Bulk concentration kg
Gc = input(nspec+2:2*nspec+1)';   %Gas concentration kg m-3
Qin = DSC.Q(find(DSC.time >= t,1 ) );  %J/s


OLc( OLc <= 0.0 | ~isreal(OLc) ) = 0.0;
Gc( Gc <= 0.0 | ~isreal(Gc) ) = 0.0;


%%Caluclate System Properties
%Recalculate Liquid Volume
X(1:nspec) = OLc(1:nspec) ./ sum(OLc);  % Mass fraction in the liquid
rhol_i = sum(X.*rhol); % Mass-weighted average density (kg/m3)
liq_vol = sum(OLc) ./ rhol_i; %m3

lconstant_volume = 1;

if lconstant_volume
    vap_vol = DSC.tot_vol - liq_vol;
    
    %Pressure at the DSC temperature
    moles_vap = DSC.moles_air + sum(Gc ./ MW .* vap_vol);
    press = moles_vap ./ vap_vol .* R .* T;
else   %Constant Pressure
    press = 101325;  %Pa
    
    %Volume at the DSC Temperature
    vap_vol = DSC.moles_air ./ (press ./R./T - sum(Gc./MW)); %m3
end


%%Calculate Intensive Properties
% Equilibrium pressures at the DSC temperature (should be Pa)
psat(1:nspec) = pstar.*exp(dHvap.*(1./T_ref - 1./T)./R);
csat(1:nspec) = MW.*psat./R./T;
% Diffusion coefficients of the species
D_air(1:nspec) = Dn_air(1:nspec).*(T/T_ref).^mu1(1:nspec);

% Calculating the composition    
OLc_tot = sum(OLc);
if sum(OLc) > 1.0e-22;
    
    % Mole fractions
    n_apu(1:nspec) = X(1:nspec) ./ MW(1:nspec);
    n_tot_apu = sum(X(1:nspec) ./ MW(1:nspec));
    Xm(1:nspec) = n_apu(1:nspec) ./ n_tot_apu;
    
    % Bulk density (Mass-weighted)
    rhol = sum(X(1:nspec).*rhol(1:nspec));
    
    % Bulk surface tension (Mole-weighted)
    sigmal = sum(Xm(1:nspec).*sigma1(1:nspec));
    
    % Bulk volume and size
    vp = OLc_tot ./ rhol;
    rp = (3.*vp./4./pi).^(1./3);
    
    % Equilibrium pressures, mole based
    peq(1:nspec) = Xm(1:nspec).*psat(1:nspec);
    
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
% pv_i(nspec) = 0;

% Mass flux of each compound between the organic liquid and gas phases   kg s-1
OL_flx = zeros(1,nspec);
for i = 1:nspec
    ln_fact = (1.0 - pv_a(i)./press)/(1.0 - pv_i(i)./press);
    if OLc_tot > 0.0 
        
        if ln_fact > 0.0  %Evaporation
                          %Stefan (Convective) Flow + Diffusive Flow
            
            OL_flx(i) = DSC.Area.*D_air(i)./DSC.L .*MW(i).*press./R./T.* ...
                log((1.0 - pv_a(i)./press)./(1.0 - pv_i(1,i)./press)); %kg/s
            
        else        %Condensation
                    %Just Diffusive flow
            
            OL_flx(i) = DSC.Area.*D_air(i)./DSC.L .* ...
                (pv_i(i)-pv_a(i)) ./R./T.*MW(i); %kg/s
        end
        
    else    %Not enough mass; Do Nothing
        OL_flx(i) = 0.0;
    end
end

if OL_flx(nspec) > 0.0
    OL_flx(nspec);
end

% Mass flux of each compound to the gas phase
g_flx(1:nspec) = -OL_flx ./ vap_vol; %kg m-3 s-1

% Calc Water Equilibrium
mass_eqh2o_vap = peq(nspec) .* vap_vol ./ R ./ T .* MW(nspec); %kg
mass_eqh2o_liq = DSC.h2o - mass_eqh2o_vap; %kg


%Energy Balance (Q is going into the system)
dT_latent = sum( OL_flx(1:nspec) .* dHvap(1:nspec) ./ MW(1:nspec) ); %J s-1
dT_Cpliq = sum( OLc(1:nspec) .* Cp_liq(1:nspec) ./ MW(1:nspec) ); %J K-1
dT_Cpvap = sum( Gc(1:nspec) .* DSC.vap_vol .* Cp_vap(1:nspec)./ MW(1:nspec)); %J K-1
dTdt = ( Qin + dT_latent) / (dT_Cpliq + dT_Cpvap);  %K s-1


%Store [Pressure, Press_H2O, Pvap_H2O, EqVapMass_H2O, EqLigMass_H2O,
%       Flux_H2O, TotLiqVol, TotVapVol]
p_h2o = [press, pv_i(nspec), pv_a(nspec), ...
    mass_eqh2o_vap, mass_eqh2o_liq, OL_flx(nspec), liq_vol, vap_vol,...
    dT_latent];
fprintf(1,'time: %4.2f  T(degC): %3.2f  H2O_flux(mg/s): %10.3e   LatentHeat: %10.3e   Lat/Q: %10.3e\n',...
    t, T-273.15, OL_flx(nspec).*1e6,dT_latent,dT_latent/Qin);



if dTdt ~= dTdt
    pause
end

% Changing to column vector
flx0 = [dTdt, OL_flx, g_flx]';
flx0( ~isreal(flx0) | isnan(flx0) ) = 0.0;




end