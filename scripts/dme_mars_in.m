% dme_mars_in.m 
%   Create input structure for drag-modulated entry at Mars. Vehicle based
%   on idea of landing an MER-class payload with MSL-class accuracy,
%   subject to current launch vehicle maximum diameter constraints
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
%   *Created AUG 2012, Z.R.Putnam

function in = dme_mars_in()
%#codegen

%% Define input structure

in = default_in;


%% Planetary data

in.p = get_planetary_data( 3, 3 ); % Load settings for Mars, exponential atmosphere


%% Simulation data

% % Termination conditions
in.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in.s.term.or.value(1) = -2000; % double, m, terminate at -2 km MOLA
in.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

% in.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
% in.s.term.or.var(1) = uint8(4); % uint8, nd, termination mode (altitude)
% in.s.term.or.value(1) = 460; % double, m, terminate at -2 km MOLA
% in.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

% Trajectory (MER-B entry conditions)
% v0 = 5500; % m/s
% fpa0 = -11.47; % deg

in.s.traj.rate = 10; % double, Hz, integration rate
in.s.traj.t_ini = 0; % double, s, initial time
in.s.traj.t_max = 2000; % double, s, maximum allowable time

% state
in.s.traj.alt = 125e3; % initial altitude, m
in.s.traj.lat = 0; % initial latitude, rad
in.s.traj.lon = 0; % initial longitude, rad
in.s.traj.vel_pp_mag = 5500; % initial velocity, m/s
in.s.traj.gamma_pp = -11*pi/180; % initial flight-path, deg->rad
in.s.traj.az = 90*pi/180; % initial azimuth, deg->rad
[in.s.traj.r_pci_ini, in.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
    in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
    in.p.omega, 0, 0, in.p.r_e, in.p.r_p);

% in.s.traj.r_pci_ini = [(in.p.r_e + 125e3); 0; 0]; % double(3), m, initial PCI position vector, 125 km altitude
% in.s.traj.v_pci_ini = v0*[sind(fpa0); cosd(fpa0); 0]; % double(3), m/s, initial PCI inertial velocity vector


%% Vehicle data

in.v.mp.m_ini = 830; % double, kg, vehicle mass at entry

% Aerodynamics
in.v.aero.mode = uint8(2); % uint8, nd, aerodynamics mode (table look-up)
in.v.aero.area_ref = pi*(4.5/2)^2; % double, m^2, aerodynamic reference area, assumes 4.7 diameter possible
in.v.aero.nose_radius = 1; % double, m, nose radius
temp = load(['.' filesep 'data' filesep 'aero_sphere_cone_70deg.mat']);
in.v.aero.table = temp.table;

% GNC
in.v.gnc.g.p.bank_mode = uint8(1); % uint8, nd, constant bank
in.v.gnc.g.p.const_bank = double(0);


end % dme_mars_in()