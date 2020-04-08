% earth_in.m 
%   basic input structure for earth entry
%   no guidance, control, etc
%
% Aeroassist Simulation source code
%
% Developed by:
%   Entry systems Design Lab (EsDL)
%   University of Colorado Boulder
%
% Inputs:
%   none
%   
% Outputs:
%   in - struct, multi, input structure with parameters for SEJ guidance
%
% Created Jan 2019, E. Roelke
% 
function in = earth_basic_in()
%#codegen

% define constants
max_alt = 125e3;
min_alt = 0;
vatm = 10e3; % m/s
efpa = -5.2; % deg
beta1 = 20;
cd = 1.05;
cl = 0;
rc = 2;
rn = rc/2;
m = beta1*cd*pi*rc^2;

%% Define input structure
in = define_in;

%% Planetary data
in.p = get_planetary_data( 2, 3 ); % Load settings for Earth with winds
temp = load('./data/atm_data/atm_earth_gram2007_nominal.mat');
in.p.atm.table = temp.table;

%% Constants
in.c.earth_g = 9.81;

%% Simulation data
% Termination conditions
in.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in.s.term.or.value(1) = min_alt; % double, m, terminate at 0m altitude
in.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

in.s.term.or.type(2) = int8(0); % int8, nd, crossing/equality condition
in.s.term.or.var(2) = uint8(1); % uint8, nd, termination mode (altitude)
in.s.term.or.value(2) = max_alt; % double, m, terminate at 125km altitude
in.s.term.or.direction(2) = int8(1); % int8, nd, positive crossing

% Trajectory
in.s.traj.rate = 100; % double, Hz, integration rate
in.s.traj.t_ini = 0; % double, s, initial time
in.s.traj.t_max = 1000; % double, s, maximum allowable time
in.s.traj.alt = max_alt;
in.s.traj.r_pci_ini = [(in.p.r_e + max_alt); 0; 0]; % double(3), m, initial PCI position vector, 125 km altitude
in.s.traj.v_pci_ini = vatm*[sind(efpa); cosd(efpa); 0]; % double(3), m/s, initial PCI inertial velocity vector
in.s.data_rate = 1;

%% Vehicle data
% Aerodynamics
in.v.aero.mode = uint8(1); % uint8, nd, aerodynamics mode
in.v.aero.aoa = 0; % double, rad, angle-of-attack
in.v.aero.cl = cl; % double, nd, lift coefficient
in.v.aero.cd = cd; % double, nd, drag coefficient
in.v.aero.area_ref = pi*rc^2; % double, m^2, aerodynamic reference area
in.v.aero.nose_radius = rn; % double, m, nose radius
in.v.mp.m_ini = m;

% Guidance
in.v.gnc.g.p.rate = 0.01; % nd, guidance rate
in.v.gnc.g.p.bank_mode = uint8(1); %  nd, constant bank
in.v.gnc.g.p.aoa_mode = uint8(1); %  nd, constant angle-of-attack
in.v.gnc.g.p.prop_mode = uint8(0); %  no guidance

% Control
in.v.gnc.c.p.rate = 0.01;
in.v.gnc.c.p.aoa_mode = uint8(0);
in.v.gnc.c.p.bank_mode = uint8(0);

% Navigation
in.v.gnc.n.p.rate = in.s.traj.rate;

end % sej_a_earth_in()
