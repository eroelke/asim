% est_atm.m
%   estimate/save atmospheric parameters for predictor guidance
%   -density corrector
%   -temp
%   -pressure
%   
% Entry systems Design Lab (EsDL)
%   University of Colorado Boulder
% 
% Written by: Evan Roelke, Dec 2018
% 	Feb 2019, add density estimation capture table, E. Roelke
% 
function guid = est_atm(in,guid,nav,calcs)
%#codegen
if (guid.p.atm.Kflag)
    atm_mode = guid.p.atm.mode; %atmospheric estimation mode
    atm_nom = guid.p.planet.atm_nom;
    rp = guid.p.planet.r_p;
    re = guid.p.planet.r_e;
    % get mass, geometry
    switch (guid.p.dm_mode)
        case 4  % dej_n
            j_ind = guid.dej_n.s.stage+1;  % current jettison stage
            
            if (j_ind > guid.dej_n.p.n_jett)
                guid.s.atm.K_dens = 1;
                return;
            end
            area_ref = guid.dej_n.p.area_refs(j_ind);  % get current stage area (1-indexed)
            mass = guid.dej_n.p.masses(j_ind);
            cd = guid.dej_n.p.cds(j_ind);
            
        case 5 % cvdma
            area_ref = guid.cvdma.s.Aref;
            cd = guid.cvdma.s.cd;
            mass = guid.cvdma.p.mass;
            
        otherwise
            guid.p.atm.Kflag = uint8(0);    % turn off density corrector
            fprintf('WARNING: No Density Corrector Input Information in Guidance Struct...\nTurning Density Corrector Off\n');
            guid.s.atm.K_dens = 1;
            return;
    end
    
    % get nav data
    V_pf_mag = norm(nav.s.v_pf_pci);    % Planet-fixed velocity (PCI), m/s
    A_mag = norm(nav.s.a_sens_pci); % Sensed acceleration (PCI), m/s^2  
    
    % get dens corrector params
    K_gain = guid.p.atm.K_gain;
    K_bounds = guid.p.atm.K_bounds;
    
    % Estimate density
    dens_est = (2*mass*A_mag) / ( (V_pf_mag^2)*area_ref*cd ); % kg/m^3
    [lat,~] = get_lat_lon(nav.s.r_pcpf,re,rp);
    alt = get_alt(nav.s.r_pcpf,re,rp,lat);
    
    % get nominal atmospheric density
    dens_model = lin_interp(atm_nom(:,1), atm_nom(:,2),alt); % kg/m^3
    
    % Density multiplier
    K_dens = dens_est/dens_model; % nd
    
    % if very large K_dens (eg. on jettison time step), keep same
    % sort of defeats the purpose of the low-pass filter but meh
    if K_dens > 2*K_bounds(2) && guid.p.dm_mode == uint8(4)
        K_dens = guid.s.atm.K_dens;
    end
    
    % Limit to min/max values
    K_dens = median([K_bounds(1), K_dens, K_bounds(2)]); % nd
    
    % low-pass filter
    guid.s.atm.K_dens = K_gain*K_dens + (1 - K_gain)*guid.s.atm.K_dens; % nd
%     guid.s.K_dens = K_dens;

    % save truth value
    guid.s.atm.rho_est = interp1(in.p.atm.table(:,1),in.p.atm.table(:,2),alt);
    guid.s.atm.K_true = guid.s.atm.rho_est / dens_model;

    % save atmospheric data
    if (atm_mode > 0)
        if alt <= guid.s.atm.atm_hist(guid.s.atm.atm_ind,1)
            % this alt hasnt saved atm data yet
            T = calcs.T + normrnd(0,in.s.sigs.T_sens);
            pres = calcs.pres + normrnd(0,in.s.sigs.P_sens);

            % rewrite alt at this step
            guid.s.atm.atm_hist(guid.s.atm.atm_ind,:) = [alt dens_est T pres];
            guid.s.atm.atm_ind = guid.s.atm.atm_ind + 1;    %index++
        else
            % on upswing out of atm
            %  -> horizontal atmospheric modeling
        end

        % get rss error for atmospheric estimates
        if (guid.p.atm.rss_flag) %ensemble filter
            guid.s.atm.rss_K = 0;
            guid.s.atm.rss_ens = 0;
            
            atm_true = flip(in.p.atm.table(:,1:2));   % true density table
            if (guid.p.atm.mode == uint8(2))
                atm_ens = flip(guid.s.atm.atm_curr);
            else
                atm_ens = nan;
            end
            h_ind = nan(1,1);
            h_ind(1) = find(atm_true(:,1) <= 70e3,1,'first');
            for i = 1:h_ind
%                 rho_true = interp1(atm_true(:,1),atm_true(:,2),alt);
%                 rho_nom = interp1(atm_nom(:,1),atm_nom(:,2),alt);
%                 rho_ens = interp1(atm_ens(:,1),atm_ens(:,2),alt);
                
                guid.s.atm.rss_K = guid.s.atm.rss_K + ... 
                    (K_dens*atm_nom(i,2) - atm_true(i,2))^2;
                
                if (guid.p.atm.mode == uint8(2))
                    if (atm_ens(1,1) == 0)
                        % ensemble not done yet, use k_dens
                        guid.s.atm.rss_ens = guid.s.atm.rss_K;
                    else
                        guid.s.atm.rss_ens = guid.s.atm.rss_ens + ...
                            (atm_ens(i,2) - atm_true(i,2))^2;
                    end
                else
                    guid.s.atm.rss_ens = 0;
                end
                
                % get nominal atm rss err if not done yet (const value)
                if guid.s.atm.rss_nom == 0
                    guid.s.atm.rss_nom = guid.s.atm.rss_nom + ...
                        (atm_nom(i,2) - atm_true(i,2))^2;
                end
            
            end %for loop  
            
            if (guid.s.atm.rss_nom ~= 0)
                guid.s.atm.rss_nom = sqrt(guid.s.atm.rss_nom)/h_ind;
            end
            guid.s.atm.rss_ens = sqrt(guid.s.atm.rss_ens)/h_ind;
            guid.s.atm.rss_K = sqrt(guid.s.atm.rss_K)/h_ind;
        else
            guid.s.atm.rss_nom = nan;
            guid.s.atm.rss_K = nan;
            guid.s.atm.rss_ens = nan;
        end %check save rss
        
    end
    
else
    guid.s.atm.K_dens = 1;
end

end % est_atm
