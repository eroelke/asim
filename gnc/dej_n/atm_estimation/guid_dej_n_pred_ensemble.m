% guid_sej_a_pred.m 
%   Predictor for single-event jettison aerocapture guidance algorithm
%   Propogates the entire trajectory forward including jettison event(s)
%   given jettison time in tj(i)
%
% Aeroassist Simulation source code
% Developed by:
%   Space Systems Design Lab, Georgia Tech
%   Colorado Center for Astrodynamics Research, CU Boulder
%
% Inputs:
%   y0 - double(6), multi, initial state [position; velocity]
%   t0 - double(1), s, initial time
%   s - NPC state structure
%   p - NPC input structure
%   j_ind - current jettison index for actual vehicle state
%   delT - small jettison time change (for newton method)
%   atm_curr - atmospheric estimate model
% 
% Outputs:
%   r_a - estimated apoapsis radius
%   dr_a - error in estimated apoapsis radius
%
% Major Revision History:
%   *Created August 2018, E. Roelke
%   *Fixed oe calcs Sep 2018 - E. Roelke
% 
function [r_a, dr_a, pflag, rss] = guid_dej_n_pred_ensemble(y0, t0, s, p, j_ind, delT, atm_curr, atm_true)
%#codegen

% Initialize
y = y0; % multi, initial state     y = [r_inertial;v_inertial]
% t = t0; % s, initial time
pflag = 0; % nd, termination flag
% h = p.t_inc_max; % s, (max) increment step size
tstep = 1/p.traj_rate;    %s, integration time step
% h = 0.1;

rss = 0;

% initialize
tj0 = s.tj_next(j_ind);        % jettison time this index
tj0 = tj0 + delT;
jett_flag = logical(false); % predictor jettison flag

i = 1;
% K_dens = s.atm.K_dens;

ttvec = t0:tstep:p.t_max;

% for debugging
i = 1;
h = nan(length(ttvec),1);
% rhoK = altitude;
% % rho_full = altitude;
% rho_model = altitude;
% rho_nom = altitude;
% rho_true = altitude;
% rhoE = altitude;
% rhoF = rhoE;
% rhoR = rhoF;
% temp = altitude;
% pres = altitude;


% time = altitude;
% time(1) = t0;

t = nan(length(ttvec)+1,1);
t(1) = t0;

% x(:,i) = y0;
% ttemp(i) = t0;
%% Integration (RK4)
while ~pflag

    
    if t(i) >= tj0 && ( jett_flag == logical(false) )
        jett_flag = logical(true);
        j_ind = j_ind + 1;
    end
        
    
    % get area, mass of current stage
    area_ref = p.area_refs(j_ind);
    mass = p.masses(j_ind);
    cd = p.cds(j_ind);
    cl = p.cls(j_ind);
    
    %{      
    testing output
    [k1,alt,rhoK(i),rho_model(i),rho_nom(i),rho_true(i)] = ... 
        npc_eom_atm(t, y, area_ref, atm_model, atm_true, K_dens, p, s, cd, cl, mass);
    altitude(i) = alt;
    i = i + 1;
    time(i) = time(i-1) + h;
    %}
    
    [k1,alt] = eom_npc_ensemble( t(i), y, area_ref, p, s, cd, cl, mass, atm_curr);
    [k2,~] = eom_npc_ensemble(t(i) + tstep/2, y + (tstep*k1/2), area_ref, p, s, cd, cl, mass, atm_curr);
    [k3,~] = eom_npc_ensemble(t(i) + tstep/2, y + (tstep*k2/2), area_ref, p, s, cd, cl, mass, atm_curr);
    [k4,~] = eom_npc_ensemble(t(i) + tstep, y + k3*tstep, area_ref, p, s, cd, cl, mass, atm_curr);
    
    y = y + (tstep/6)*(k1 + 2*k2 + 2*k3 + k4); % set new state vector
    
    if (i == length(t))
        pflag = 3;
        break;
    end
    t(i+1) = t(i) + tstep;
        
%     rho_true(i) = rho_dat.rho_true;
%     rhoK(i) = rho_dat.rho_K;
%     rho_model(i) = rho_dat.rho_model;
%     rho_nom(i) = rho_dat.rho_nom;
%     rhoE(i) = rho_dat.rhoE;
%     rhoF(i) = rho_dat.rhoF;
%     rhoR(i) = rho_dat.rhoR;
%     rho_full(i) = rho_dat.rho_full;
    
    h(i) = alt;

    %% Check termination conditions
    if alt > p.planet.alt_max
        pflag = 1; % Maximum altitude exceeded
        break;
    elseif alt < p.planet.alt_min
        pflag = 2; % Minimum altitude exceeded
        break;
%     elseif t > p.t_max
%         pflag = 3; % Maximum time exceeded
%         break;
    end
    
    i = i + 1;
end % while ~plag

% calc orbit elements
r = y(1:3); % m, position vector
v = y(4:6); % m/s, inertial velocity vector
oe = rv2oe([r;v],p.planet.mu);    % get orbital elements from r, v
a = oe(1);
e = oe(2);

% apoapsis estimate
r_a = a*(1+e);       % apoapsis for given jettison time
dr_a = r_a - (p.tgt_ap + p.planet.r_e);   %apoapsis radius error, m

% 
% if (t0 > 50)
%     keyboard;
% end

% keyboard;
%     rho_nom = p.planet.atm_table(:,2);
%     RvarT = nan(1000,1);
%     RvarM = RvarT;
% %     RvarK = RvarT;
%     for j = 1:length(atm_true)
% %         RvarK(j) = rhoK(j)/rho_nom(j);      % density corrector variation
%         RvarM(j) = atm_curr(j,2)/rho_nom(j); % model variation
%         RvarT(j) = atm_true(j,2)/rho_nom(j);  % actual variation
% %         RvarF(j) = rhoF(j)/rho_nom(j);
% %         RvarR(j) = rhoR(j)/rho_nom(j);
%         
% %         RvarFull(j) = rho_full(j)/rho_nom(j);
%         
% %         rhoExp = 
% %         RvarExp(j) = 
% %         RvarR(j) = rho_ref(j)/rho_nom(j);
%         
% %         diff(j) = rho_model(j) - RvarM(j);
% %         rss = rss + (diff(j)^2);
% %         diff(j) = K_dens(j) - Rvar(j);
% %         rss = rss + ( diff(j)^2 );
%           altitude(j) = atm_true(j,1);
%     end
% 
%     figure(); hold on
%     grid on
%     set(gca,'FontSize',14)
%     title(['Predictor Ensemble Density Variation vs. Altitude, Tol=' num2str(p.atm.tol)],'FontSize',10)
%     plot(RvarM,altitude./1000,'b','LineWidth',2)
% %     plot(RvarS,altitude(1:min_ind)./1000,'g','LineWidth',2)
% %     plot(RvarR(1:min_ind),altitude(1:min_ind)./1000,'g','LineWidth',2)
%     plot(s.atm.K_dens.*ones(length(RvarT),1),altitude./1000,'r','LineWidth',2)
%     plot(RvarT,altitude./1000,'k--','LineWidth',2)
% %     plot(RvarFull(1:min_ind),altitude(1:min_ind)./1000,'c','LineWidth',2)
%     plot(s.atm.K_dens,(norm(y0(1:3))-p.planet.r_e)/1000,'*r','LineWidth',2)
%     xlabel('Density Variation (kg/m^3)')
%     ylabel('Altitude (km)');
%     ylim([min(altitude)/1000 - 3 max(altitude)/1000 + 1])
%     legend('Correlation Filter','Density Factor','Truth', ... 
%         'Predictor Initiation','location','se','FontSize',10)



end % guid_sej_a_pred


