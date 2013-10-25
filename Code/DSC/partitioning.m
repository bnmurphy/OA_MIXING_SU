function [X_eq Xm_eq caer_tot_eq rp_eq ctot] = ...
    partitioning(caer_tot, delta_ROG, fi,ctot_i, rp, sigma,MW, rho, csat, nspec, temp, ntot)
 
    if isempty(ctot_i)
     % Total mass of each compound kg/m3
     ctot = fi.*delta_ROG./0.08206./temp.*(MW.*1e3)./1e9;   
    else
     ctot = ctot_i;
    end
     % Total moles of each compound
     ctot_moles = MW.*ctot;
     % Initial guess for mass fractions
     X = ones(1,nspec)./nspec;
     R = 8.314472;
     
     tol = 1.0;
     caer_tot_new = 1;
     
     if ~isempty(rp)
         
         % The option when the radius of the particle is known (initial equilibration)
         
        while(tol>1e-9 & caer_tot_new > 0)
        % Total mass of each species in the aerosol
        mspec = X.*caer_tot_new;
        % Mole fraction of each species in the aerosol
        Xm = (X./MW)./sum(X./MW);
        % Density
        rho_ave = sum(X.*rho);
        % Surface tension
        sigma_ave = sum(Xm.*sigma);
        % Kelvin effect
        Ke = exp(2.*sigma_ave.*MW./rho_ave./rp./R/temp);
        % Equilibrium concentration. NOTE! These are done with mass
        % fractions in the original yield fits!
        %ceq = Xm.*csat.*Ke;
            ceq = X.*csat.*Ke;
        % Fraction of each compound in the aerosol phase
        for i = 1:nspec
            if mspec(i) > 0
            Xp(i) = 1./(1 + ceq(i)./mspec(i));
            else
                Xp(i) = 1.0;
            end
        end
        % Mass of each compound in the aerosol phase
        caer = Xp.*ctot;
        % Total aerosol concentration
        caer_tot_new = sum(caer);
       if caer_tot_new > 0.0
%        New mass fractions
        X = caer./caer_tot_new;
        tol = abs((caer_tot_new - caer_tot))./caer_tot;
        caer_tot = caer_tot_new;
        else
 %       X = caer.*0.0;
        tol = 0.0;
        caer_tot_new = 0.0;
        end
        caer_tot = caer_tot_new;
        end
        
     elseif ~isempty(ntot)
         
         % The option when the total number concentration of the particles
         % is known (final equilibration)
         
         while(tol>1e-6 & caer_tot_new > 0)
        % Mass of each individual particle
        mp = caer_tot_new./ntot;
        % Total mass of each species in the aerosol
        mspec = X.*caer_tot_new;
        % Mole fraction of each species in the aerosol
        Xm = (X./MW)./sum(X./MW);
        % Density
        rho_ave = sum(X.*rho);
        % Volume of each particle
        vp = mp./rho_ave;
        % Radius of each particle
        rp = (3.*vp/4./pi).^(1./3);
        % Surface tension
        sigma_ave = sum(Xm.*sigma);
        % Kelvin effect
        Ke = exp(2.*sigma_ave.*MW./rho_ave./rp./R/temp);
        % Equilibrium concentration. NOTE! These are done with mass
        % fractions in the original yield fits!
        ceq = X.*csat.*Ke;
                % Fraction of each compound in the aerosol phase
         for i = 1:nspec
            if mspec(i) > 0
            Xp(i) = 1./(1 + ceq(i)./mspec(i));
            else
                Xp(i) = 1.0;
            end
         end
        % Mass of each compound in the aerosol phase
        caer = Xp.*ctot;
        % Total aerosol concentration
        caer_tot_new = sum(caer);
       if caer_tot_new > 0.0
%        New mass fractions
        X = caer./caer_tot_new;
        tol = abs((caer_tot_new - caer_tot))./caer_tot;
        caer_tot = caer_tot_new;
        else
        X = caer.*0.0;
        tol = 0.0;
        caer_tot_new = 0.0;
        end
        end 
        caer_tot = caer_tot_new;
     end
        
        X_eq = X;
        Xm_eq = (X_eq./MW)./sum(X_eq./MW);
        caer_tot_eq = caer_tot;
        mp_eq = caer_tot_eq./ntot;
        if mp_eq > 0.0
        rho_eq = sum(X_eq.*rho);
        else
        rho_eq = 1000;
        end
        vp_eq = mp_eq./rho_eq;
        rp_eq = (3.*vp_eq/4./pi).^(1./3);
        
end