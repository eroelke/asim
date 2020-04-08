% Major Revision History:
%   Earlier - Functionality written, Zach Putnam
%   July 2017 - Script developed, Michael Werner
%   February 2018 - Formulated into wrapper script, Evan Roelke
% 
% function to use bisection method to calculate the required FPA corridor
% for a given aeroshell geometry and target apoapsis
% Inputs:
%   planet - string of planet eg. 'venus'
%   x_atm: entry state vector
%       (1) - h_atm, km
%       (2) - v_atm, m/s
%       (3) - fpa_atm (nominal), deg
%       (4) - lat_atm, deg
%       (5) - lon_atm, deg
%       (6) - az_atm, deg
%   Aeroshell aerodynamics
%       cd: 1x2 vector of pre and post jettison cd vals
%       cl: 1x2 vector of pre and post jettison cd vals
%   Geometry:
%       m: 1x2 vector of pre and post jettison masses, kg
%       rn: 1x2 vector of pre and post jettison nose radii, m
%       rc: 1x2 vector of pre and post jettison backshell radii, m
%   Targeting:
%       ha_tgt: target apoapsis altitude, m
%       ha_tol: tolerance on target apoapsis altitude, m
%       fpa_tol: tolerance on final FPA corridor, deg
%       corr0: initial guess on fpa corridor - must bracket the solution

function [beta1,beta2,beta_ratio,fpa_corr] = ... 
    calc_fpa_corr(planet,x0,cd,cl,m,rn,rc,ha_tgt,ha_tol,fpa_tol,corr0)


%% Setup - Basic Inputs
% Vehicle properties
% R = 0.75;           % Pre-jettison radius (m)
% r = 0.175;          % Post-jettison radius (m)
% area = [pi*(R)^2, pi*(r)^2];    % Drag area (m^2)
% m1 = 80;            % Pre-jettison mass (kg)
% m2 = 40;            % Post-jettison mass (kg)
% rn = 1/2*r;         % Nose radius (m)
% delta = 45;         % Sphere-cone half angle (deg)
% cd = 1.05;

%% Setup - Input Structure
in0 = get_planet_struct(planet);
% in0 = dme_venus_in; % nominal input structure

% Initial state: Trajectory properties
in0.s.traj.gamma_pp = x0.efpa*pi/180; %entry flight path angle (nominal), rad
in0.s.traj.alt = x0.h0;      % entry altitude, m
in0.s.traj.lat = x0.lat0*pi/180;      % initial latitude, rad
in0.s.traj.lon = x0.lon0*pi/180;      % initial longitude, rad
in0.s.traj.az = x0.az0*pi/180;       % initial azimuth, deg->rad
% in0.s.traj.vel_pp_mag = x_atm(2);   %entry velocity, m/s
in0.s.traj.vel_pp_mag = find_vel_pp_mag(x0.v0, in0.s.traj.gamma_pp, in0.p.omega, ... 
    (in0.p.r_e + in0.s.traj.alt)); % initial velocity, m/s

% Convert scalars to inertial position and velocity vectors
[in0.s.traj.r_pci_ini, in0.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in0.s.traj.lat, in0.s.traj.lon, in0.s.traj.alt, ...
    in0.s.traj.gamma_pp, in0.s.traj.az, in0.s.traj.vel_pp_mag, ...
    in0.p.omega, in0.s.traj.t_ini, 0, in0.p.r_e, in0.p.r_p);

% Vehicle properties
area = pi.*[rc(1)^2 rc(2)^2];       % drag areas pre and post jettison
in0.v.mp.m_ini = m(1);              % Initial mass (kg)
in0.v.mp.m_jettison = m(1)-m(2);    % Jettisoned mass (kg)
in0.v.aero.area_ref = area(1);      % Area before and after jettison (m^2)
in0.v.aero.mode = uint8(1);         % User-specified lift and drag coefficients
in0.v.aero.cd = cd(1);              % pre jettison drag coeff
in0.v.aero.cl = cl(1);              % pre jettison lift coeff
in0.v.aero.nose_radius = rn(1);     % Initial nose radius (m)

% Sim settings
in0.s.traj.rate = 50;       % double, Hz, integration rate
in0.s.traj.t_ini = 0;       % double, s, initial time
in0.s.traj.t_max = 2000;    % double, s, maximum allowable time

% Guidance Settings
min_alt = 60e3;            % minimum allowable altitude (m)
max_alt = x0.h0;     % maximum allowable altitude (m) - same as h_atm

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
% fpa_guess = [-15 -5]*pi/180;    % rad, entry FPA guess for bisection (MUST bracket corridor)
% tol_fpa = 0.01*pi/180; % rad, tolerance on final FPA corridor
% tol_ha = 1e4; % m, tolerance on final apoapse altitude

%% Calculations - Find Corridor
% Case 1: Pre-jettison, minimum beta (shallow side of corridor)
in(1) = in0;

% Case 2: Post-jettison, maximum beta (steep side of corridor)
in(2) = in0;
in(2).v.aero.area_ref = area(2);
in(2).v.mp.m_ini = m(2);
in(2).v.aero.cd = cd(2);
in(2).v.aero.cl = cl(2);
in(2).v.aero.nose_radius = rn(2);

% Compute ballistic coefficients and ratio of post/pre jettison
beta1 = in(1).v.mp.m_ini/(in(1).v.aero.area_ref*in(1).v.aero.cd);   %pre jettison
beta2 = (in(2).v.mp.m_ini)/(in(2).v.aero.area_ref*in(2).v.aero.cd); %post jettison
beta_ratio = beta2/beta1;

% Compute corridor
% fpa_tol = fpa_tol*pi/180;   %deg->rad
% corr0 = corr0.*pi/180;
[fpa_corr(1), ~, ~, ~, ~] = find_fpa_given_ha(corr0, ha_tgt, in(1), fpa_tol, ha_tol); % shallow side
[fpa_corr(2), ~, ~, ~, ~] = find_fpa_given_ha(corr0, ha_tgt, in(2), fpa_tol, ha_tol); % steep side

% fpa_corr = fpa_corr.*180/pi;    %convert corridors back to degrees

% Output corridor
% fprintf('Shallow side: %4.2f deg,\n', fpa(1)*180/pi);
% fprintf('Steep side: %4.2f deg.\n', fpa(2)*180/pi);
% 
% fprintf('Corridor width: %4.2f deg.\n', abs(fpa(2)-fpa(1))*180/pi);


end




function in = get_planet_struct(planet)

if ischar(planet) == true
    models = {'venus','earth','mars','jupiter','saturn','titan','uranus','neptune'};
    planet = find(strcmpi(planet,models) == 1);
end

switch planet
    case 1  %mercury
        error('No Planet Input Struct for Mercury\n')
    case 2  %venus
        in = dme_venus_in;
    case 3  %earth
        in = dme_earth_in;
    case 4  %mars
        in = dme_mars_in;
    case 5  %jupiter
        error('No Planet Input Struct for Jupiter\n')
    case 6  %saturn - titan
        fprintf('Planet Input Structure for Titan\n')
        in = dme_titan_in;
    case 7  %uranus
        error('No Planet Input Struct for Uranus\n')
    case 8  %neptune
        error('No Planet Input Struct for Neptune\n')
    
    otherwise
        error('Incorrect Planet Input')
end

end
