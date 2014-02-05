function [flx0, p_h2o] = fluxes_bulk(t,input,dt,nspec,Cv_sld,Cp_liq,Cv_vap,...
    pstar,dHfus,dHvap,Tmelt,T_ref,MW,sigma1,rhol,D_liq,Dn_air,mu1,alpha_m,...
    DSC,molec_group_flag,molec_group_stoich)

R = 8.314472;
%Load Local Variables
T = input(1);  %Kelvin
Sc = input(2:nspec+1)';    %Solid phase bulk mass [kg]
OLc = input(nspec+2:2*nspec+1)';  %Organic Liquid Bulk mass [kg]
Gc = input(2*nspec+2:3*nspec+1)';   %Gas concentration [kg m-3]
Qin = DSC.Q(find(DSC.time >= t,1 ) );  %J/s

Sc( Sc <= 0.0 | ~isreal(Sc) ) = 0.0;
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
Sc_tot = sum(Sc);
if sum(OLc) > 1.0e-22;
    
    % Mole fractions
    n_apu(1:nspec) = X(1:nspec) ./ MW(1:nspec);
    n_tot_apu = sum(X(1:nspec) ./ MW(1:nspec));
    Xm(1:nspec) = n_apu(1:nspec) ./ n_tot_apu;
    
    %Get Activity Coefficients
    ln_gamma = AIOMFAC_gamma_SR_v1(nspec, Xm, molec_group_flag, molec_group_stoich, T);
    gamma = exp(nonzeros(ln_gamma))';  %Convert activity coeff sparse matrix to vector    
%     gamma = ones(1,nspec);
    
    % Bulk density (Mass-weighted)
    rhol = sum(X(1:nspec).*rhol(1:nspec));
    
    % Bulk surface tension (Mole-weighted)
    sigmal = sum(Xm(1:nspec).*sigma1(1:nspec));
    
    % Bulk volume and size
    vp = OLc_tot ./ rhol;
    rp = (3.*vp./4./pi).^(1./3);
    
    % Equilibrium pressures, mole based
    peq(1:nspec) = gamma(1:nspec) .* Xm(1:nspec) .* psat(1:nspec);
    
else
    mp = 0.0;
    vp = 0.0;
    rp = 0.0;
    peq(1:nspec) = zeros(1,nspec);
    beeta(1:nspec) = zeros(1,nspec);
end

% Pressure at the bulk surface (Pa)
pv_a(1:nspec) = peq;

% Partial pressures far away from the bulk surface (Pa)
pv_i(1:nspec) = Gc(1:nspec).*R.*T./MW(1:nspec);
% pv_i(nspec) = 0;


%%%%%%%%%%%%%%%%%%%
% Gas-Liquid Mass flux  kg s-1
%%%%%%%%%%%%%%%%%%%
G2OL_flx = zeros(1,nspec);
for i = 1:nspec
    ln_fact = (1.0 - pv_a(i)./press)/(1.0 - pv_i(i)./press);
    if OLc_tot > 0.0 
        
        if ln_fact > 0.0  %Calculate Flux to the Particles            
            G2OL_flx(i) = DSC.Area.*D_air(i)./DSC.L .*MW(i).*press./R./T.* ...
                log((1.0 - pv_a(i)./press)./(1.0 - pv_i(1,i)./press)); %kg/s
        else
            G2OL_flx(i) = DSC.Area.*D_air(i)./DSC.L .* ...
                (pv_i(i)-pv_a(i)) ./R./T.*MW(i); %kg/s
        end
        
    else    %Not enough mass; Do Nothing
        G2OL_flx(i) = 0.0;
    end
end

if G2OL_flx(nspec) > 0.0
    G2OL_flx(nspec);
end

%%%%%%%%%%%%%%%%%%%
% Solid-Liquid Mass flux  kg s-1
%%%%%%%%%%%%%%%%%%%
%Calcuate solubility (xsol) as a fxn of T
xsol(1:nspec) = 1./gamma(1:nspec) .* exp(dHfus./R .* (1./Tmelt - 1./T));

Sc2OL_flx = zeros(1,nspec);
% for i = 1:nspec
%     ln_fact = (1.0 - xsol)/(1.0 - Xm);
%     if Sc_tot > 0.0 
%         
%         if ln_fact > 0.0  %Calculate flux to the solid phase
%             Sc2OL_flx(i) = DSC.Area.*D_liq(i)./(DSC.L./10) .*MW(i).*rhol.* ...
%                 log((1.0 - xsol(i))./(1.0 - Xm(i))); %kg/s
%             
%         else        
%             Sc2OL_flx(i) = DSC.Area.*D_liq(i)./(DSC.L./10) .*MW(i).*rhol.* ...
%                 (Xm(i)-xsol(i)); %kg/s
%         end
%         
%     else    %Not enough mass; Do Nothing
%         Sc2OL_flx(i) = 0.0;
%     end
% end
% 
% if Sc2OL_flx(nspec) > 0.0
%     Sc2OL_flx(nspec);
% end



% Sum up all of the production and loss terms
Sc_flx(1:nspec) = +Sc2OL_flx;  %kg s-1
OLc_flx(1:nspec) = G2OL_flx - Sc2OL_flx; %kg s-1
Gc_flx(1:nspec) = -G2OL_flx ./ vap_vol; %kg m-3 s-1

% Calc Water Equilibrium
mass_eqh2o_vap = peq(nspec) .* vap_vol ./ R ./ T .* MW(nspec); %kg
mass_eqh2o_liq = DSC.h2o - mass_eqh2o_vap; %kg


%Energy Balance (Q is going into the system)
dT_latent_gasliq = sum( G2OL_flx(1:nspec) .* dHvap(1:nspec) ./ MW(1:nspec) ); %J s-1
dT_latent_solliq = 0.0; %sum( Sc2OL_flx(1:nspec) .* dHfus(1:nspec) ./ MW(1:nspec) ); %J s-1

% dT_Cvsol = sum( Sc(1:nspec) .* Cv_sld(1:nspec) ./ MW(1:nspec) ); %J K-1
dT_Cpliq = sum( OLc(1:nspec) .* Cp_liq(1:nspec) ./ MW(1:nspec) ); %J K-1
dT_Cvvap = sum( Gc(1:nspec) .* DSC.vap_vol .* Cv_vap(1:nspec)./ MW(1:nspec)); %J K-1
dTdt = ( Qin + dT_latent_gasliq) / ... % + dT_latent_solliq) / ...
        (dT_Cpliq + dT_Cvvap); % + dT_Cvsol);  %K s-1


    
%Store [Pressure, Press_H2O, Pvap_H2O, EqVapMass_H2O, EqLigMass_H2O,
%       Flux_H2O, TotLiqVol, TotVapVol]
p_h2o = [press, pv_i(nspec), pv_a(nspec), ...
    mass_eqh2o_vap, mass_eqh2o_liq, G2OL_flx(nspec), liq_vol, vap_vol,...
    dT_latent_gasliq, dT_latent_solliq];
% fprintf(1,'time: %4.2f  T(degC): %3.2f  H2O_gas_liq_flux(mg/s): %10.3e   LatentHeat: %10.3e   Lat/Q: %10.3e\n',...
%     t, T-273.15, G2OL_flx(nspec).*1e6,dT_latent_gasliq,dT_latent_gasliq/Qin);
fprintf(1,'time: %4.2f  T(degC): %3.2f  OLc_flux(mg/s): %10.3e   LatentHeat: %10.3e   dTCpliq: %10.3e   dTCpVap: %10.3e   dTdt: %10.3e\n',...
    t, T-273.15, G2OL_flx(1).*1e6,dT_latent_gasliq,dT_Cpliq,dT_Cvvap,dTdt);


if dTdt ~= dTdt
    pause
end

% Changing to column vector
flx0 = [dTdt, Sc_flx, OLc_flx, Gc_flx]';
flx0( ~isreal(flx0) | isnan(flx0) ) = 0.0;




end