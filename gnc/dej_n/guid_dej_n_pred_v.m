% guid_dej_n_pred_v.m 
%   Predictor for discrete-event jettison. Saves state throughout predictor
%   for debugging / plotting / etc.
% 
% Warning: slow!
%
% Aeroassist Simulation source code
% Entry systems Design Lab (EsDL)
% University of Colorado Boulder
%
% Inputs:
%   y0 - double(6), multi, initial state [position; velocity]
%   t0 - double(1), s, initial time
%   s - NPC state structure
%   p - NPC input structure
%   j_ind - current jettison index for actual vehicle state
%   dt - small jettison time change (for newton method)
%
% Outputs:
%   r_a - estimated apoapsis radius
%   dr_a - error in estimated apoapsis radius
%
% *Created Sep 2018, E. Roelke
% *Added outside loop for multiple stages, Mar 2019, E. Roelke    
% 
function [r_a,dr_ap,dr_ap_true,pflag] = guid_dej_n_pred_v(y0, t0, s, p, j_ind, dt, guid)
%#codegen

% Initialize
% y = y0; % multi, initial state     y = [r_inertial;v_inertial]
% t = t0; % s, initial time
% h = p.t_inc_max; % s, (max) increment step size
h = 1/p.traj_rate;    %s, integration time step

K_dens = guid.s.atm.K_dens;

% for debugging
% i = 1;
% x(:,i) = y0;
% ttemp(i) = t0;

y = y0;
t = t0;
jett_flag = logical(false);
j_ind0 = j_ind;
pflag = 0;
tj0 = s.tj_next(j_ind);        % jettison time this index
tj0 = tj0 + dt;

% atm_model = guid.s.atm_est;
% atm_true = in.p.atm.table
%initialize some vars for the compiler
rho_true = nan; 
rho_di = rho_true; rho_K = rho_true; rho_nom = rho_true;
K_nom = rho_true;
K_model = K_nom;
model_err = K_nom;
model_var = K_nom;
K_truth = K_nom;
true_var = K_nom;
rho_ecf = K_nom;
rho_hybrid = K_nom;

y(:,1) = y0;
alt = norm(y0(1:3)) - guid.p.planet.r_e;
vmag = norm(y0(4:6));
t = t0;
idx = 1;
%% Integration (RK4)
while ~pflag
    
    if t(idx) >= tj0 && ( jett_flag == logical(false) )
        jett_flag = logical(true);
        j_ind = j_ind + 1;
    end
        
    % get area, mass of current stage
    area_ref = p.area_refs(j_ind);
    mass = p.masses(j_ind);
    cd = p.cds(j_ind);
    cl = p.cls(j_ind);
    
    % integrate with rk4
    [k1,alt(idx),rho_dat] = npc_eom_atm(t(idx), y(:,idx), area_ref, K_dens, p, s, cd, cl, mass, guid);
    [k2,~,~] = npc_eom_atm(t(idx) + h/2, y(:,idx) + (h*k1/2), area_ref, K_dens, p, s, cd, cl, mass, guid);
    [k3,~,~] = npc_eom_atm(t(idx) + h/2, y(:,idx) + (h*k2/2), area_ref, K_dens, p, s, cd, cl, mass, guid);
    [k4,~,~] = npc_eom_atm(t(idx) + h, y(:,idx) + k3*h, area_ref, K_dens, p, s, cd, cl, mass, guid);
    y(:,idx + 1) = y(:,idx) + (h/6)*(k1 + 2*k2 + 2*k3 + k4); % set new state vector
    
    rho_true(idx) = rho_dat.rho_true;   % actual density
    rho_di(idx) = rho_dat.rho_di; % model density (interpolator or ecf)
    rho_K(idx) = rho_dat.rho_K; % scale factor
    rho_nom(idx) = rho_dat.rho_nom; % nominal density
    rho_ecf(idx) = rho_dat.rho_ecf; % ecf density
    rho_hybrid(idx) = rho_dat.rho_hybrid;   % hybrid model density
    
    %% Check termination conditions
    if alt(idx) > guid.p.planet.alt_max
        pflag = 1; % Maximum altitude exceeded
        break;
    elseif alt(idx) < guid.p.planet.alt_min
        pflag = 2; % Minimum altitude exceeded
        break;
%     elseif t > p.t_max
%         pflag = 3; % Maximum time exceeded
%         break;
    elseif (isnan(alt(idx)))
        pflag = 3;
        break;
    end
    
    t(idx + 1) = t(idx) + h;
    vmag(idx + 1) = norm(y(4:6,idx));
    idx = idx + 1;
end % while ~plag


% if (guid.p.atm.mc_ind == guid.s.atm.ind_curr - 1)
% keyboard;
% 
% 
% model = rho_ecf;
% 
% K_nom = rho_nom / rho_true;
% for i = 1:length(model)
%     K_model(i) = model(i) / rho_nom(i);
%     model_var(i) = (model(i) - rho_nom(i)) / rho_nom(i);
%     K_truth(i) = lin_interp(guid.p.planet.atm_true(:,1),guid.p.planet.atm_true(:,2),alt(i)) / rho_nom(i);
%     true_var(i) = (rho_true(i) - rho_nom(i)) / rho_nom(i);
% %     model_err(i) = (rho_di(i) - rho_true(i)) / rho_true(i);
% end
% 
% figure(); hold on
% grid on;
% set(gca,'FontSize',14)
% title(['Predictor Initiated At: ' num2str(round(alt(1)/1000,0)) ' km']);
% plot(K_model, alt./1000,'LineWidth',2);
% plot(K_truth, alt./1000,'k--','LineWidth',2)
% xline(1,'Color',[0 0 0] + 0.4,'LineWidth',1.5)
% plot(K_truth(1), alt(1)./1000,'b*','LineWidth',1.5)
% yline(alt(1)/1000)
% xlabel('Density Variation');
% ylabel('Altitude (km)');
% legend('Hybrid','Truth','location','ne')
% 
% 
% % figure(); hold on
% % semilogy(guid.p.planet.atm_true(:,2)); hold on
% % semilogy(guid.s.atm.atm_curr(:,2))
% 
% end

r = y(1:3,end); % m, position vector
v = y(4:6,end); % m/s, inertial velocity vector
[r_a, ~, e] = get_apoapsis_radius([r;v], guid.p.planet.mu); % apoapsis for given jettison time

dr_ap_true = r_a - (p.tgt_ap + guid.p.planet.r_e);   %apoapsis radius error, m
dr_ap = dr_ap_true - p.bias(j_ind0).tgt_ap; %minus because negative error is below target

if (e > 1)
    pflag = 3;  % hyperbolic traj
    return;
end

end % guid_sej_a_pred
