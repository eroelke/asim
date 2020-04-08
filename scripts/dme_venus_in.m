% dme_venus_in.m 
%   Create input structure for drag-modulated entry at Venus.
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   none
%   
% Outputs:
%   in - struct, multi, input structure with default values
%
% Major Revision History:
%   *Created OCT 2015, M.S.Werner

function in = dme_venus_in()
%#codegen

%% Define input structure

in = define_in;


%% Planetary data

in.p = get_planetary_data( 1, 3 ); % Load settings for Venus, GRAM atmosphere + winds


%% Simulation data

% % Termination conditions
in.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in.s.term.or.value(1) = 0; % double, m, terminate at 0m altitude
in.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

in.s.traj.rate = 10; % double, Hz, integration rate
in.s.traj.t_ini = 0; % double, s, initial time
in.s.traj.t_max = 2000; % double, s, maximum allowable time

% state
in.s.traj.alt = 125e3; % initial altitude, m
in.s.traj.lat = 0; % initial latitude, rad
in.s.traj.lon = 0; % initial longitude, rad
in.s.traj.vel_pp_mag = 5500; % initial velocity, m/s
in.s.traj.gamma_pp = -5*pi/180; % initial flight-path, deg->rad
in.s.traj.az = 90*pi/180; % initial azimuth, deg->rad
[in.s.traj.r_pci_ini, in.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
    in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
    in.p.omega, 0, 0, in.p.r_e, in.p.r_p);


%% Vehicle data
% Nominal Values based on 'Aerocapture CubeSat Demo Study Report'

in.v.mp.m_ini = 25; % double, kg, vehicle mass at entry

% Aerodynamics
in.v.aero.mode = uint8(1); % uint8, nd, aerodynamics mode (user-provided)
in.v.aero.area_ref = pi*(1.2/2)^2; % double, m^2, aerodynamic reference area
in.v.aero.cl = 0;
in.v.aero.nose_radius = 0.025; % double, m, nose radius
temp = load(['./data/aero_data/aero_sphere_cone_45deg.mat']);
in.v.aero.table = temp.table;

% GNC
in.v.gnc.g.p.bank_mode = uint8(1); % uint8, nd, constant bank
in.v.gnc.g.p.const_bank = double(0);


end