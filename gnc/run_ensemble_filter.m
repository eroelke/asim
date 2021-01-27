% run_ensemble_filter.m
%   select atmospheric structure from all disperse atmospheres
% 
% Entry systems Design Laboratory (EsDL)
% University of Colorado Boulder
% 
% Written By: Evan Roelke, Jul 2019
% 
% Inputs:
%   p: guidance parameter input struct
%   s: guidance state struct
% 
% atm_all.mat is structured as follows:
%   col 1: altitude
%   cols 2-1001: GRAM density profiles
%   
% thus, any GRAM index that is found is 2-indexed
% 
function [s] = run_ensemble_filter(y0, t0, p, s, dej_n)
%#codegen
if (p.atm.Kflag && p.atm.mode >= uint8(3))
    % speed up process
%     if (s.atm.ind_curr-1 ~= p.atm.mc_ind)
        model_ind = int32(find(isnan(s.atm.atm_hist(:,2)) == true, 1, 'first') - 1);
        model_ind = model_ind(1);   % for compiler to know it has dimension 1x1
        rho_est = s.atm.atm_hist(model_ind,2);

        tol = rho_est*p.atm.ens_tol;
        
        % C compiler does not like loading dynamic
        % string so this is the way to do it...
        switch (p.planet.p_ind)
            case 1 %venus
                atm = load('./data/atm_data/atm_venus_all.mat');
            case 2 %earth
                atm = load('./data/atm_data/atm_earth_all.mat');
            otherwise
                fprintf('Warning: No suitable atmospheric table for ensemble correlation filter.\n');
                return;
        end
        atms = atm.mcs;
        alt = alt_from_guid(y0, t0, p.planet);
        
        h_ind = find(atms(:,1) >= alt,1,'first');   % get alt from nominal profile
        mcs = atms(h_ind,2:end)';    % all MC density values this altitude
        
        % awkward coding to allow it to compile
        temp = [];
        while (isempty(temp))
            upper = rho_est + tol;
            lower = rho_est - tol;
            % get indices of gram profiles with densities
            % within search tolerance at this alt
            inds = find(mcs >= lower & mcs <= upper);   
            if ~isempty(inds)
                break;
            end
            tol = tol + p.atm.ens_tol;
        end
        
        rss_rho = zeros(length(inds),1);        
        rss_true = rss_rho;
        
        for k = 1:length(inds)
            atm_curr = [atms(:,1) atms(:,inds(k)+1)];    % monte carlo profile for this index
            for jj = 1:model_ind
                % density rss
                rho_curr = lin_interp(atm_curr(:,1), atm_curr(:,2), s.atm.atm_hist(jj,1)); %interp to get mc atm density at atm model alt
                diff = (s.atm.atm_hist(jj,2) - rho_curr)^2;
                rss_rho(k) = rss_rho(k) + diff;
                
                rho_true = interp1(p.planet.atm_true(:,1), p.planet.atm_true(:,2), s.atm.atm_hist(jj,1));
                diff_true = (rho_curr - rho_true)^2;
                rss_true(k) = rss_true(k) + diff_true;
                
                % pressure rss
                % temp rss
            end
        end % for each atm mc profile
        
        
        %%% make a graphic detailing search tolerance
        %{
        atm_model = s.atm.atm_hist;
        x_min = (alt/1000) - 0.5;
        x_max = (alt/1000) + 0.5;
        y_min = log(lower)-0.05;
        y_max = log(upper)+0.15;
        figure(); hold on
        set(gca,'FontSize',14)
        for j = 1:1000
            if mod(j,10) == 0  %only look at some
                if isempty(find(j == inds,1))
                    % not in set
                    p8=plot(atms(h_ind,1)/1000,log(atms(h_ind,j+1)),'b.','MarkerSize',15);
                    if mod(j,randi(75)) == 0
                        plot(atms(:,1)./1000,log(atms(:,j+1)),'b--','LineWidth',1.5)
                    end
                else
                    p7=plot(atms(h_ind,1)/1000,log(atms(h_ind,j+1)),'r.','MarkerSize',15);
                    if mod(j,randi(75)) == 0
                        plot(atms(:,1)./1000,log(atms(:,j+1)),'r--','LineWidth',1.5)
                    end
                end
            end
        end
        p1=plot(atms(h_ind,1)/1000,log(rho_est),'k.','MarkerSize',25);
        p2=plot(atms(h_ind,1)/1000,log(upper),'.','Color',[0 0 0]+0.5,'MarkerSize',20);
        p3=plot(atms(h_ind,1)/1000,log(lower),'.','Color',[0 0 0]+0.5,'MarkerSize',20);
        p4=plot(atm_model(:,1)/1000,log(atm_model(:,2)),'k--','LineWidth',2);
        p5=plot(p.planet.atm_true(:,1)/1000,log(p.planet.atm_true(:,2)),'--','Color',[0 0.5 0],'LineWidth',1.5);
        p6=plot(p.planet.atm_true(h_ind,1)/1000,log(p.planet.atm_true(h_ind,2)),'.','Color',[0 0.5 0],'MarkerSize',15);
        xlabel('Altitude (km)')
        ylabel('Log \rho (Log kg/m^3)')
        xlim([x_min x_max])
        ylim([y_min y_max])
        set( gca, 'xdir', 'reverse' );
        legend([p1,p2,p3,p7,p8,p4,p5], ...
            'Current \rho_{est}', ...
            'Upper \rho Search Bound', ...
            'Lower \rho Search Bound', ...
            '\rho Profile Within Bounds', ...
            '\rho Profile Outside Search', ...
            'Density History, \rho_{hist}', ...
            'True Density Profile, \rho_{true}','Location','eastoutside');
        
        %}
        
        % get lowest rss error (density)
        s.atm.ind_curr = nan(1,1);
        if (t0 < (dej_n.tj_curr(dej_n.stage + 1) - 5))
            s.atm.ind_curr(1) = inds(rss_rho == min(rss_rho)) + 1;
        else
            s.atm.ind_curr(1) = find(s.atm.ecf_scores == ... 
                max(s.atm.ecf_scores),1);
        end
        s.atm.ecf_scores(s.atm.ind_curr) = s.atm.ecf_scores(s.atm.ind_curr) + 1;
        
%         s.atm.atm_curr = [atms(:,1) s.atm.K_dens .*atms(:,s.atm.ind_curr)];
        s.atm.atm_curr = [atms(:,1) atms(:,s.atm.ind_curr)];
%         ind_true = find(rss_true == 0);
%         sig_curr = (atms(h_ind,s.atm.ind_curr) - rho_est)/atms(h_ind,s.atm.ind_curr);
        
%     else
%         s.atm.atm_curr(:,2) = s.atm.K_dens .* p.planet.atm_true(:,2);
%     end % speed up
    
else
    s.atm.atm_curr = p.planet.atm_nom(:,1:2);
end
end %run_ensemble_filter