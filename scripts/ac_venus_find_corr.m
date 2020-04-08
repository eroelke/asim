% ac_venus_find_corr.m
%   function to determine corridor width of aerocapture trajectory at earth
% 
% Entry systems Design Laboratory
% University of Colorado Boulder
% 
% Written By:  Evan Roelke
% 
% Inputs:
%   aero:       aerodynamics and geometry input structure
%       m:          mass vector (kg) [2x1]
%       rc:         backshell radii (m) [2x1]
%       cd:         drag coefficient [1x1] - const across jettison
%   v_atm:      planet-relative entry speed (km/s)
%   ha_tgt:     target apoapsis altitude (km)
% 
% Outputs:
%   fpa_corr:   entry corridor width (deg)

function fpa_corr = ac_venus_find_corr(aero,v_atm,ha_tgt)

%% Setup - Basic Inputs
% Vehicle properties
R = aero.rc(1);           % Pre-jettison radius (m)
r = aero.rc(2);          % Post-jettison radius (m)
area = [pi*(R)^2, pi*(r)^2];    % Drag area (m^2)
m1 = aero.m(1);            % Pre-jettison mass (kg)
m2 = aero.m(2);            % Post-jettison mass (kg)
rn = 1/2*r;         % Nose radius (m)
delta = 45;         % Sphere-cone half angle (deg)
cd = aero.cd;

%% Setup - Input Structure
in0 = dme_venus_in; % nominal input structure
% temp = load('C:\Users\Evan\Documents\Academics\Research\asim\data\atm_venus_mc.mat');
% temp = load('C:\Users\Evan\Documents\Academics\Research\asim\data\atm_venus_JPL.mat');
[nom, ~] = parse_atm_data(in0, 'venus');
in0.p.atm.table = nom;

% Initial state: Trajectory properties
in0.s.traj.gamma_pp = -5.5*pi/180;
in0.s.traj.alt = 150e3; % initial altitude, m
in0.s.traj.lat = 0*pi/180; % initial latitude, rad
in0.s.traj.lon = 0*pi/180; % initial longitude, rad
in0.s.traj.az = 90*pi/180; % initial azimuth, deg->rad
in0.s.traj.vel_pp_mag = v_atm*1e3;

% Convert scalars to inertial position and velocity vectors
[in0.s.traj.r_pci_ini, in0.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in0.s.traj.lat, in0.s.traj.lon, in0.s.traj.alt, ...
    in0.s.traj.gamma_pp, in0.s.traj.az, in0.s.traj.vel_pp_mag, ...
    in0.p.omega, in0.s.traj.t_ini, 0, in0.p.r_e, in0.p.r_p);

% Vehicle properties
in0.v.mp.m_ini = m1; % Initial mass (kg)
in0.v.mp.m_jettison = m1-m2; % Jettisoned mass (kg)
in0.v.aero.area_ref = area(1); % Area before and after jettison (m^2)
in0.v.aero.mode = uint8(1);     % User-specified lift and drag coefficients
in0.v.aero.cd = cd;
in0.v.aero.cl = 0;
in0.v.aero.nose_radius = rn; % Nose radius (m)

% Sim settings
in0.s.traj.rate = 200; % double, Hz, integration rate
in0.s.traj.t_ini = 0; % double, s, initial time
in0.s.traj.t_max = 2000; % double, s, maximum allowable time

% Guidance Settings
min_alt = 0;        % minimum allowable altitude (m)
max_alt = 150e3;    % maximum allowable altitude (m)
in0.v.gnc.g.p.dm_mode = uint8(0);

% Termination conditions
in0.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(1) = min_alt; % double, m, terminate at 0m altitude
in0.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

in0.s.term.or.type(2) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(2) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(2) = max_alt; % double, m, terminate at 125km altitude
in0.s.term.or.direction(2) = int8(1); % int8, nd, positive crossing

%% Setup - User Inputs

% ha_tgt = 500e3;    % m, target post-aerocapture apogee altitude
% ha_tgt = 400e3;

fpa_guess = [-11 -0.1];    % deg, entry FPA guess for bisection (MUST bracket corridor)

tol_fpa = 0.01; % deg, tolerance on final FPA corridor
tol_ha = 25; % km, tolerance on final apoapse altitude

%% Calculations - Find Corridor
% Case 1: Pre-jettison, minimum beta (shallow side of corridor)
in(1) = in0;
in(1).v.aero.area_ref = area(1);

% Case 2: Post-jettison, maximum beta (steep side of corridor)
in(2) = in0;
in(2).v.aero.area_ref = area(2);
in(2).v.mp.m_ini = m2;

% Compute ballistic coefficients
beta_1 = in(1).v.mp.m_ini/(in(1).v.aero.area_ref*cd);
beta_2 = (in(2).v.mp.m_ini)/(in(2).v.aero.area_ref*cd);

% Compute corridor
[fpa(1), ~, ~, ~, ~] = find_fpa_given_ha(fpa_guess, ha_tgt, in(1), tol_fpa, tol_ha);
[fpa(2), ~, ~, ~, ~] = find_fpa_given_ha(fpa_guess, ha_tgt, in(2), tol_fpa, tol_ha);

% [fpa(1)] = fpa_given_haf(-5.2, ha_tgt, in(1), tol_fpa, tol_ha);
% [fpa(2)] = fpa_given_haf(-5.2, ha_tgt, in(2), tol_fpa, tol_ha);

fpa_corr = abs(fpa(1)-fpa(2)); %deg

% Output corridor
% fprintf('For Target Apoapsis %4.1f km:\n',ha_tgt/1000)
% fprintf('Shallow side: %4.4f deg,\n', fpa(1)*180/pi);
% fprintf('Steep side: %4.4f deg.\n', fpa(2)*180/pi);
% 
fprintf('Corridor width: %4.4f deg.\n', fpa_corr);

end
