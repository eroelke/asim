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
%   atm_model: density estimates from density estimator func (est_K_dens.m)
%   j_ind - current jettison index for actual vehicle state
%   dt - small jettison time change (for newton method)
%
% Outputs:
%   r_a - estimated apoapsis radius
%   dr_a - error in estimated apoapsis radius
%
% Major Revision History:
%   *Created August 2018, E. Roelke
%   *Fixed oe calcs Sep 2018 - E. Roelke
% 
% E. Roelke - I don't recommend using this at all. It is mainly for
% personal research and is a nightmare of code attempts and ideas.
% 
function [r_a,dr_a,pflag,atm_err] = guid_dej_n_pred_atm(y0, t0, s, p, atm_model, atm_true, j_ind, dt, guid)
%#codegen

% Initialize
y = y0; % multi, initial state     y = [r_inertial;v_inertial]
t = t0; % s, initial time
pflag = 0; % nd, termination flag
% h = p.t_inc_max; % s, (max) increment step size
h = 1/p.traj_rate;    %s, integration time step
% h = 0.1;

tj0 = s.tj_next(j_ind);        % jettison time this index
tj0 = tj0 + dt;

% prevent prediction at previous time
% while tj0 <= t0
%     tj0 = tj0 + 5;
% end

jett_flag = logical(false); % predictor jettison flag

i = 1;
K_dens = s.atm.K_dens;

ttvec = t0:h:p.t_max;

% for debugging
i = 1;
altitude = nan(length(ttvec),1);
rhoK = altitude;
% rho_full = altitude;
rho_model = altitude;
rho_nom = altitude;
rho_true = altitude;
rhoE = altitude;
rhoF = rhoE;
rhoR = rhoF;
temp = altitude;
pres = altitude;


% time = altitude;
% time(1) = t0;

t = nan(length(ttvec)+1,1);
t(1) = t0;

% x(:,i) = y0;
% ttemp(i) = t0;
%% Integration (RK4)
while ~pflag

    % get jettison estimate
    
    % only project the current stage and next jettison event instead of
    % all?
    
    % save delta_ap at each jettison, and if other stage cannot improve upon last stage by
    % the time it should jettison, do not jettison
        % logic for improving on n-2+ jettison stages?
    
%     if j_ind <= p.n_jett
%         tj0 = s.tj_next(j_ind);        % jettison time this index
%     else
%         tj0 = s.tj_next(j_ind-1);
%     end
 

    % Set step size to ensure step at jettison
%     if (tj0 >= t) && (tj0 <= (t + h) )
%         h = tj0 - t;  % inner time step
%         if h == 0
%             h = p.t_inc_min;
%         end
%         if t >= tj0
% %             jett_flags(j_ind) = true;
%             j_ind = j_ind + 1;          % next jettison stage
%         end
%     else
%         h = p.t_inc_max;    % step size
% %         h = p.t_inc_min;
%     end
    
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
    
    [k1,alt,rho_dat] = npc_eom_atm(t(i), y, area_ref, atm_model, atm_true, ... 
        K_dens, p, s, cd, cl, mass, guid);
    [k2,~,~] = npc_eom_atm(t(i) + h/2, y + (h*k1/2), area_ref, atm_model, ... 
        atm_true, K_dens, p, s, cd, cl, mass, guid);
    [k3,~,~] = npc_eom_atm(t(i) + h/2, y + (h*k2/2), area_ref, atm_model, ... 
        atm_true, K_dens, p, s, cd, cl, mass, guid);
    [k4,~,~] = npc_eom_atm(t(i) + h, y + k3*h, area_ref, atm_model, ... 
        atm_true, K_dens, p, s, cd, cl, mass, guid);
    
    y = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4); % set new state vector
    
    if (i == length(t))
        pflag = 3;
        break;
    end
    t(i+1) = t(i) + h;
    
    rho_true(i) = rho_dat.rho_true;
    rhoK(i) = rho_dat.rho_K;
    rho_model(i) = rho_dat.rho_model;
    rho_nom(i) = rho_dat.rho_nom;
%     rhoE(i) = rho_dat.rhoE;
%     rhoF(i) = rho_dat.rhoF;
%     rhoR(i) = rho_dat.rhoR;
%     rho_full(i) = rho_dat.rho_full;
    
    altitude(i) = alt;
    
    temp(i) = rho_dat.T;
    pres(i) = rho_dat.P;
        
%     if (p.atm_mode == uint8(2))
%         K_dens(i+1) = K_new;
%         K_true(i) = rho_dat.rho_true/rho_dat.rho_nom;
%     else
%         K_dens(i+1) = K_dens(1);
%     end
        
    % get altitude
%     alt = norm(y(1:3)) - p.planet.r_e;
    
    % save state (for debugging)
%     i = i + 1;
%     x(:,i) = y;  % ensure jettison event(s) test
%     ttemp(i) = t;
%     m(i) = mass;  

    %% Check expedience conditions
%     if alt <= p.alt_min
%         h = 2; % going to impact, expediate process
%     end

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

% check to see if jettison events worked - for debugging
% r = x(1:3,:);
% v = x(4:6,:);
% for i = 1:size(x,2)
%     vmag(i) = norm(v(:,i));
%     rmag(i) = norm(r(:,i));
% end
% plot(vmag./1000,rmag./1000)
% plot(ttemp,vmag./1000)

% if p.atm_mode ~= uint8(0)
%     figure();hold on
% %     plot(rhoK,altitude./1000,'r');
% %     plot(rho_model,altitude./1000,'b');
%     set(gca,'FontSize',14)
% %     plot(rho_model - rhoK,altitude./1000,'LineWidth',2)
%     plot(
%     xlabel('\rho_{model} - \rho_{K_{dens}}')
%     ylabel('Altitude (km)')
%     legend(['Starting Altitude: ' num2str(altitude(1)/1000) ' km'],'location','nw')
% end
% calc orbit elements and params - get orbital apoapsis error from inertial
% state
r = y(1:3); % m, position vector
v = y(4:6); % m/s, inertial velocity vector
oe = rv2oe([r;v],p.planet.mu);    % get orbital elements from r, v
a = oe(1);
e = oe(2);

r_a = a*(1+e);       % apoapsis for given jettison time
dr_a = r_a - (p.tgt_ap + p.planet.r_e);   %apoapsis radius error, m



% figure(); hold on
% plot(rho_true,'r')
% plot(rho_est,'b')
% plot(rho_nom,'k')

%     if t0 > 80
%         keyboard
%     end
%     

%     idx = find(isnan(rho_ref) == 1,1)-1;

%     figure();
%     grid on
%     set(gca,'FontSize',14)
%     semilogx(rho_ref(1:idx),altitude(1:idx)./1000,'r','LineWidth',1.5); hold on
%     semilogx(rho_model,altitude./1000,'c','LineWidth',1.5); hold on
%     semilogx(rhoK,altitude./1000,'b','LineWidth',1.5); hold on
%     semilogx(rho_true,altitude./1000,'--k','LineWidth',1.5); hold off
%     set(gca,'FontSize',14)
%     legend('Reference Trajectory','Atmospheric Capture','Density Corrector','Truth','location','ne')
%     xlabel('Density (kg/m^3)')
%     ylabel('Altitude (km)')

%  if (t0 > 60)
%      keyboard
%  end
 
%  figure(); hold on
%  plot(rho_model)
%  plot(rho_true)
%  plot(rhoK)
%  
%  min_ind = find(isnan(altitude) == true,1,'last');

% visualize RSS error between K_dens and model
%     min_ind = find(altitude == min(altitude));
% %     min_ind = find(isnan(altitude)==true,1)-1;
%     RvarK = nan(min_ind,1);
%     RvarM = RvarK;
%     RvarT = RvarK;
% %     RvarR = RvarK;
% %     diff = RvarM;
% %     RvarF = RvarK;
%     RvarR = RvarK;
% %     RvarFull = RvarK;
%     rss = 0;
%     for j = 1:min_ind
%         RvarK(j) = rhoK(j)/rho_nom(j);      % density corrector variation
%         RvarM(j) = rho_model(j)/rho_nom(j); % model variation
%         RvarT(j) = rho_true(j)/rho_nom(j);  % actual variation
% %         RvarF(j) = rhoF(j)/rho_nom(j);
%         RvarR(j) = rhoR(j)/rho_nom(j);
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
%     end
%     
%     figure(); hold on
%     grid on
%     set(gca,'FontSize',14)
%     title('Predictor Density Variation vs. Altitude, Earth','FontSize',10)
%     plot(RvarM,altitude(1:min_ind)./1000,'b','LineWidth',2)
% %     plot(RvarS,altitude(1:min_ind)./1000,'g','LineWidth',2)
%     plot(RvarR(1:min_ind),altitude(1:min_ind)./1000,'g','LineWidth',2)
%     plot(RvarK(1:min_ind),altitude(1:min_ind)./1000,'r','LineWidth',2)
%     plot(RvarT(1:min_ind),altitude(1:min_ind)./1000,'k','LineWidth',2)
% %     plot(RvarFull(1:min_ind),altitude(1:min_ind)./1000,'c','LineWidth',2)
%     plot(RvarK(1),altitude(1)/1000,'*r','LineWidth',2)
%     xlabel('Density Variation (kg/m^3)')
%     ylabel('Altitude (km)');
%     ylim([min(altitude)/1000 - 3 max(altitude)/1000 + 1])
%     legend('Atmospheric Capture','Density Corrector','Truth', ... 
%         'Predictor Initiation','location','ne','FontSize',10)

% figure(); hold on
% % plot(t(1:length(altitude)),temp)
% semilogy(t(1:length(altitude)),rhoK); hold on
% semilogy(t(1:length(altitude)),rho_nom); hold on
% semilogy(t(1:length(altitude)),rho_true); hold on
% semilogy(t(1:length(altitude)),rho_model); hold on

% figure(); hold on
% yyaxis left
% plot(t(1:length(altitude)),temp)
% yyaxis right
% semilogy(t(1:length(altitude)),rho_model)

% figure(); hold on
% semilogy(t,rhoF); hold on
% semilogy(t,rho_true,'k'); hold on
% semilogy(t,rho_nom);



% calculate rho_model RSS error this time step
if ~isempty(find(isnan(rho_model)==false,1))
    sig = 0;
    for i = 1:length(altitude)
        if (isnan(rho_model(i)) == false)
            sig = sig + (rho_model(i) - rho_true(i))^2;
        end
    end
    if (isnan(sig) == true)
        sig = 0;
    end
    atm_err = sqrt(sig/(length(altitude)-1));
else
    atm_err = nan;
end

end % guid_sej_a_pred



% function dx = eoms(x0)
%     r = x0(1);  %pos
%     v = x0(2);  %vel
%     theta = x0(3);  %lon
%     phi = x0(4);    %lat
%     gamma = x0(5);  %fpa
%     psi = x0(6);    %heading angle
%     p = x0(7);  %pres
%     rho = x0(8);    %dens
%     
%     
%     
%     dr = v*sin(gamma);
%     dv = (fT/m) - g*sin(gamma) + w^2*r*cos(phi) * ... 
%         ( sin(gamma)*cos(phi) - cos(gamma)*sin(phi)*sin(psi) );
%     dTheta = (v*cos(gamma)*cos(psi))/(r*cos(phi));
%     dPhi = (v*cos(gamma)*sin(psi))/r;
%     dGamma = (1/v)*
%     dPsi = (1/v)*( (fN*sin(nu))/(m*cos(gamma)) - ... 
%         ((v^2)/r)*cos(gamma)*cos(psi)*tan(phi) + ...
%         2*w*v*(tan(gamma)*cos(phi)*sin(psi) - sin(phi)) - ...
%         (w^2*r/cos(gamma))*sin(phi)*cos(phi)*cos(psi) );
%     dp = -rho0*g*v*sin(gamma);
%     dRho = -p0^2*g*v*sin(gamma)/rho0;
%     
%     xdot = [dr;dv;dTheta;dPhi;dGamma;dPsi;dp;dRho];
%     
% end


% function dy = 

