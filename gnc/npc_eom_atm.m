% npc_eom.m
%   numerical predictor-corrector force calcs
%   updated to include improved atm model
%       verbose output for figures
% 
% Entry systems Design Lab (EsDL)
% University of Colorado Boulder
% 
% Written by: Evan Roelke, Feb 2019
% 
% E. Roelke - I don't recommend using this at all. It is mainly for
% personal research and is a nightmare of code attempts and ideas.
%
function [ydot,alt,rho_dat] = ... 
    npc_eom_atm(t, y, area_ref, atm_model, ...
    atm_true, K_dens, p, s, cd, cl, mass, guid)
%#codegen

% Assign states
r_ii = y(1:3); % Inertial position
v_ii = y(4:6); % Inertial velocity
rmag_ii = norm(r_ii); % Inertial position magnitude

% planet params
omega = p.planet.omega;     % planet rotation rate (rad/s)
rp = p.planet.r_p;
re = p.planet.r_e;
mu = p.planet.mu;
J2 = p.planet.j2;

% PCI -> PCPF
theta = norm(omega)*t; % Rotation angle [rad]
L_PI = get_LPI(theta);  % DCM for PCI->PCPF

% Inertial and planet centric velocities
r_pp = L_PI*r_ii; % Position vector planet/planet [m]
v_pp = L_PI*(v_ii - cross(omega,r_ii)); % Velocity vector planet/planet [m/s]

pos_pp_mag = norm(r_pp);
vel_pp_mag = norm(v_pp);

% Angular Momentum Calculations
h_pp = cross(r_pp,v_pp);
hmag_pp = norm(h_pp);
hhat_pp = h_pp/hmag_pp;

% Compute latitude and longitude
[lat,lon] = get_lat_lon(r_pp,re,rp);
% Compute NED basis unit vectors
[uN,uE,uD] = get_unit_NED( lat, lon ); % nd

% Inertial flight path angle
% arg = median([-1, 1, h_ii_mag/(pos_ii_mag*vel_ii_mag)]); % limit to [-1,1]
% gamma_ii = acos(arg);
% if dot(pos_ii,vel_ii) < 0
%     gamma_ii = -gamma_ii;
% end

% Relative flight path angle
arg = median([-1, 1, hmag_pp/(pos_pp_mag*vel_pp_mag)]); % limit to [-1,1]
gamma_pp = acos(arg);
if dot(r_pp,v_pp) < 0
    gamma_pp = -gamma_pp;
end

% Compute Altitude
if re == rp
    alt = rmag_ii - re;
else
    alt = get_alt(r_pp,re,rp,lat);
end

% Get expected (nominal) atmospheric params
% [rho,~,~,wind] = get_atm_data(in,alt);

atm = interp1(p.planet.atm_table(:,1), p.planet.atm_table(:,2:7), alt);
rho_nom = atm(1);
rho_K = rho_nom*K_dens; % density corrector on nominal atm model (for monte carlo)    
wind = [atm(4) atm(5) atm(6)];  % east; north; vertical


% rho_mc = nan(1000,1);
% for i = 1:1000
%     rho_mc(i) = lin_interp(mc(i).table(:,1),mc(i).table(:,2),alt);
% end


% interpolate altitude against saved model to get density
idend = find(isnan(atm_model(:,2))==true,1)-1;


% curr = atm_model(1:idend,:);
% cov0 = cov(curr);   %covariance on current atmos model
% rho0 = rho_nom;   % expected density

% keyboard;

if p.atm.atm_mode > uint8(0)
    if alt >= atm_model(idend,1)

        ind = find(alt >= atm_model(:,1),1); % find first index for altitude
        ind = ind(1);   % force scalar (dumbass mex compiler...)
        if ~isempty(ind)
            if ind == 1
                rho_model = atm_model(1,2);
                T = atm_model(1,3);
                P = atm_model(1,4);
            else
                rho_model = interp1(atm_model(ind-1:ind,1),atm_model(ind-1:ind,2),alt);
                T = interp1(atm_model(ind-1:ind,1),atm_model(ind-1:ind,3),alt);
                P = interp1(atm_model(ind-1:ind,1),atm_model(ind-1:ind,4),alt);
            end
        else
            rho_model = nan;
            T = nan;
            P = nan;
        end

        % save density estimate
        if isnan(rho_model) == false
            rho = rho_model;
        else
            rho = rho_K;
        end
    else
        % kalman filter - to do

        % use density corrector
        rho = rho_K;
        rho_model = nan;
        T = nan;
        P = nan;
    end
else
    rho = rho_K;
    rho_model = nan;
    T = nan;
    P = nan;
end

% save density data (visualization purposes)
rho_true = interp1(atm_true(:,1),atm_true(:,2),alt);
rho_dat.rho_true = rho_true(1);
rho_dat.rho_model = rho_model;
rho_dat.rho_K = rho_K;
rho_dat.rho_nom = rho_nom;
rho_dat.T = T;
rho_dat.P = P;
% rho_dat = nan;


% if (p.atm.atm_mode == uint8(5))
   
%     keyboard;
    
    % spline
%     model_ind = guid.s.atm_ind-1;
%     a = spline(atm_model(1:model_ind,1),atm_model(1:model_ind,2));
%     
%     s = spline(guid.dej_n.p.planet.atm_table(:,1),guid.dej_n.p.planet.atm_table(:,2));
%     figure(); hold on
%     plot(atm_model(1:model_ind,1),log(ppval(s,atm_model(1:model_ind))))
%     plot(atm_model(1:model_ind,1),log(ppval(a,atm_model(1:model_ind))))
%     plot(atm_true(:,1),log(atm_true(:,2)));
%     plot(flip(guid.dej_n.p.planet.atm_table(:,1)),log(flip(guid.dej_n.p.planet.atm_table(:,2))),'--');
%         
%     figure(); hold on
%     plot(atm_model(:,1),log(atm_model(:,2)))
%     plot(atm_true(:,1),log(atm_true(:,2)))
    
    
    % nonlinear regression - not working
%     func = @(b,x) b(1) + b(2)*x(:,1).^b(3) + b(4)*x(:,2).^b(5);
%     beta0 = 0.2*ones(1,5);
%     model = fitnlm(atm_model(1:idend,1),atm_model(1:idend,2),func,beta0);
% end



% if (p.atm.atm_mode == uint8(3))
%     
%     % justus1969
%     [L1, L3] = getAtmL(alt,'rho');
%     L = L1*L3/sqrt( L1^2 * sin(gamma_pp)^2 + L3^2 * cos(gamma_pp)^2 );
%     
%     rho_dat.rhoE = rho_nom*(exp(-rho_nom/L));
%     
%     val = sin(2*pi*r_pp(1)/(1000*L1))*sin(2*pi*r_pp(3)/(1000*L3));
%     vari = (.004^2 / 2)*(1 - cos(2*pi*r_pp(1)/1000/L1)*cos(2*pi*r_pp(3)/1000/L3));
%     rhoF = rho_nom*(val+1);
%         
%     rho_dat.rhoF = median([rhoF - vari, rhoF, rhoF + vari]);
%     
%     
%     % justus1974
%     model = atm_model(1:idend,2);
%     R = nan(length(model)-1,1);
%     for i = 1:length(model)-1
%         R(i) = corrcoef(model(i)/rho_nom,model(i+1)/rho_nom);
% 
%         h = atm_model(i,1);
%         nom(i) = interp1(p.planet.atm_table(:,1),p.planet.atm_table(:,2),h);
%         
%         covariance = cov(atm_model(i,:),atm_model(i+1,:));
% 
%         sigRR(i) = covariance(1,1);
%         K(i) = atm_model(i,2)/nom(i);  %density var        
%     end
%     
%     R = cov(K);
%     A = R*K_dens/K(end);
%     B = K_dens*sqrt(1-R^2);
%     
%     rho_dat.rhoR = rho_nom*(A*K(end) + B);   
%     
% else
%     rho_dat.rhoE = nan;
%     rho_dat.rhoR = nan;
% end


% Convert wind to pp (PCPF) frame
wE = wind(1); % positive to the east, m/s
wN = wind(2); % positive to the north, m/s
wU = wind(3); % positive up, m/s
wind_pp = wN*uN + wE*uE - wU*uD; % v_wind/planet (PCPF), m/s
% vel_pp_rw = v_pp + wind_pp; % relative wind vector, m/s (wrong sign?)
vel_pp_rw = v_pp - wind_pp; % v_sc/wind = v_sc/p - v_wind/p
vel_pp_rw_mag = norm(vel_pp_rw); % relative wind magnitude, m/s
vel_pp_rw_hat = vel_pp_rw / vel_pp_rw_mag; % relative wind unit vector, nd

% Dynamic pressure
q = 1/2*rho*vel_pp_rw_mag^2; % base on wind-relative velocity

% Gravity force: modified Oct 2018, E. Roelke
Fg_inv_ii = -mu*mass*r_ii/(rmag_ii^3); % inverse square gravity
Fg_j2_ii = -( (3*J2*mu*re^2)/(2*rmag_ii^5) ) * ... 
    [(5*(r_ii(3)^2)/(rmag_ii^2) - 1)*r_ii(1); ...    % (5z^2/r^2 - 1)*x
     (5*(r_ii(3)^2)/(rmag_ii^2) - 1)*r_ii(2); ...    % (5z^2/r^2 - 1)*y
     (5*(r_ii(3)^2)/(rmag_ii^2) - 3)*r_ii(3) ];  % j2 gravitational accel (PCI frame)
Fg_ii = Fg_inv_ii + Fg_j2_ii;   % total PCI gravity force

% Vectors and Rotation Tensors of Interest
n1 = skew(-hhat_pp);
R1 = eye(3) + sin(pi/2)*n1 + (1-cos(pi/2))*n1*n1;
n2 = skew(vel_pp_rw_hat);
R2 = eye(3) + sin(s.bank)*n2 + (1-cos(s.bank))*n2*n2;

% Aerodynamic Forces (constant cl/cd for now)
% [cl,cd,veh] = get_aero_data(mach, veh, alt, vel_pp_rw_mag, q, cg, in); % Drag coefficient
drag_pp_hat = -vel_pp_rw_hat; % Planet relative drag force direction
drag_pp = q*cd*area_ref*drag_pp_hat; % Planet relative drag force vector
lift_pp_hat = R2*R1*vel_pp_rw_hat; % Planet relative lift force direction
lift_pp = q*cl*area_ref*lift_pp_hat; % Planet relative lift force vector
drag_ii = L_PI'*drag_pp; % Inertial drag force vector
lift_ii = L_PI'*lift_pp; % Inertial lift force vector

% Total force vector
force_ii = drag_ii + lift_ii + Fg_ii;   % no thrust or decelerator

% total acceleration on body
a_ii = force_ii/mass;

% if (p.atm.atm_mode == uint8(2))
% %     A_mag = norm(drag_ii);  % norm of drag accel
% % %     V_pf_mag = norm(v_pp);   % norm of vel
% %     dens_est = (2*mass*A_mag) / ( (vel_pp_rw_mag^2)*area_ref*cd );
% %     K_new = dens_est/rho_nom;
% %     K_new = median([guid.K_bounds(1), K_new, guid.K_bounds(2)]); % nd
% %     K_new = guid.K_gain*K_new + (1 - guid.K_gain)*K_new; % nd
% 
%     % get dynamics
%     dv = -(rho*vel_pp_rw_mag^2)/(2*beta) + g*sin(gamma)
% 
% end

% eom output
ydot = [v_ii;a_ii];
end % eom
