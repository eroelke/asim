% gnc_eom.m
%   numerical predictor-corrector force calcs
%	mirrors main integration traj_calcs to reduce predictor bias
% 
% Entry systems Design Lab (EsDL)
% University of Colorado Boulder
% 
% Written by: Evan Roelke, Oct 2018
%   Oct 2019: atmospheric estimation stuff combined
% 
function [ydot,alt,rho_dat] = npc_eom_atm(t, y, area_ref, ... 
    K_dens, p, s, cd, cl, mass, guid)
%#codegen
% Assign states
r_ii = y(1:3); % Inertial position
v_ii = y(4:6); % Inertial velocity
rmag_ii = norm(r_ii); % Inertial position magnitude

% planet params
omega = guid.p.planet.omega;     % planet rotation rate (rad/s)
rp = guid.p.planet.r_p;
re = guid.p.planet.r_e;
mu = guid.p.planet.mu;
J2 = guid.p.planet.j2;

% PCI -> PCPF
theta = norm(omega)*t; % Rotation angle [rad]
L_PI = get_LPI(theta);  % DCM for PCI->PCPF

% Inertial and planet centric velocities
r_pp = L_PI*r_ii; % Position vector planet/planet [m]
v_pp = L_PI*(v_ii - cross(omega,r_ii)); % Velocity vector planet/planet [m/s]


% Angular Momentum Calculations
h_pp = cross(r_pp,v_pp);
hmag_pp = norm(h_pp);
hhat_pp = h_pp/hmag_pp;

% Compute latitude and longitude
[lat,lon] = get_lat_lon(r_pp,re,rp);
% Compute NED basis unit vectors
[uN,uE,uD] = get_unit_NED( lat, lon ); % nd

% Compute Altitude
if re == rp
    alt = rmag_ii - re;
else
    alt = get_alt(r_pp,re,rp,lat);
end

% if (isnan(alt))
%     keyboard
% end

% Get atmospheric density params
atm = lin_interp(guid.p.planet.atm_nom(:,1), ... 
    guid.p.planet.atm_nom(:,2:7), alt); %better than interp1, handles OOB cases
rho_nom = atm(1);
wind = [atm(4) atm(5) atm(6)];  % east; north; vertical
atm_curr = guid.s.atm.atm_curr; %current atmospheric model

rho_K = rho_nom*K_dens; % density corrector on nominal atm model (for monte carlo)    
rho_true = rho_nom * guid.s.atm.K_true;
rho_di = calc_rho_interp( alt, guid.s.atm.atm_hist, K_dens*atm(1) );
rho_ecf = lin_interp( atm_curr(:,1), atm_curr(:,2), alt );

% hybrid di/ecf
atm_model = guid.s.atm.atm_hist;
% find first non-NaN density value
idend = find(isnan(atm_model(:,2)) == true, 1) - 1;
if (alt >= atm_model(idend,1))
    rho_hybrid = calc_rho_interp( alt, guid.s.atm.atm_hist, K_dens*atm(1) );
else
    atm_curr = guid.s.atm.atm_curr; %current atmospheric model
    rho_hybrid = lin_interp( atm_curr(:,1), atm_curr(:,2), alt );
end

switch (guid.p.planet.mode)
    case {0, 1}  %exponential
        wind = [0 0 0];
        rho = guid.p.planet.rho0 * exp(-alt / guid.p.planet.H);
    case 3 %table look up w winds
        atm = lin_interp(guid.p.planet.atm_nom(:,1), guid.p.planet.atm_nom(:,2:7), alt);
        rho = atm(1);   %default, expected value
        
        % get atmospheric density
        if (guid.p.atm.Kflag)
            switch (guid.p.atm.mode)
                case 1 %density factor
                    rho = rho_K;    
                case 2 %density interp
                    rho = rho_di;
                case 3 %ensemble filter
                    rho = rho_ecf;
                case 4 %DI/ECF hybrid
                    rho = rho_hybrid;
                otherwise %density factor
                    rho = rho_nom; % density corrector on nominal atm model (for monte carlo)    
            end
        end
    otherwise %exponential
        wind = [0 0 0];
        rho = guid.p.planet.rho0 * exp(-alt / guid.p.planet.H);
end

% save density data
rho_dat.rho_true = rho_true;
rho_dat.rho_di = rho_di;
rho_dat.rho_K = rho_K;
rho_dat.rho_nom = rho_nom;
rho_dat.rho_ecf = rho_ecf;
rho_dat.rho_hybrid = rho_hybrid;

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

% eom output
ydot = [v_ii;a_ii];
end % gnc_eom


