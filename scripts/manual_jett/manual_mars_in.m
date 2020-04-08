% manual_mars_in.m 
%   Create input structure for single-event jettison aerocapture at
%   Mars using manual, user-defined jettison times
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
%   in - struct, multi, input structure with parameters for SEJ DCF guidance
%
% Major Revision History:
%   *Created Jan 2017, M. Werner

function in = manual_mars_in()
%#codegen


%% Define input structure
in = define_in;


%% Planetary data

in.p = get_planetary_data( 3, 3 ); % Load settings for Mars, table atm model

%% Constants

in.c.earth_g = 9.81;


%% Simulation data

% Termination conditions
in.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in.s.term.or.value(1) = 150e3; % double, m, terminate at atmospheric exit
in.s.term.or.direction(1) = int8(1); % int8, nd, positive crossing

in.s.term.or.type(2) = int8(-1); % int8, nd, crossing/equality condition
in.s.term.or.var(2) = uint8(1); % uint8, nd, termination mode (altitude)
in.s.term.or.value(2) = 0.0; % double, m, terminate at atmospheric exit

% Trajectory
v0 = 10e3; % m/s
fpa0 = -5.2; % deg
alt0 = 150000; % m

in.s.traj.rate = 100; % double, Hz, integration rate
in.s.traj.t_ini = 0; % double, s, initial time
in.s.traj.t_max = 1000; % double, s, maximum allowable time
in.s.traj.r_pci_ini = [(in.p.r_e + alt0); 0; 0]; % double(3), m, initial PCI position vector, 125 km altitude
in.s.traj.v_pci_ini = v0*[sind(fpa0); cosd(fpa0); 0]; % double(3), m/s, initial PCI inertial velocity vector

in.s.data_rate = 1;

%% Vehicle data

beta1 = 10;
beta2 = 100;

% Aerodynamics
in.v.aero.mode = uint8(2); % uint8, nd, aerodynamics mode
in.v.aero.aoa = 0; % double, rad, angle-of-attack
in.v.aero.cl = 0.0; % double, nd, lift coefficient
in.v.aero.cd = 1.53; % double, nd, drag coefficient
in.v.aero.area_ref = 200; % double, m^2, aerodynamic reference area
in.v.aero.nose_radius = 1; % double, m, nose radius
temp = load(['.' filesep 'data' filesep 'aero_sphere_cone_70deg.mat']);
in.v.aero.table = temp.table;

in.v.mp.m_ini = beta1*in.v.aero.cd*in.v.aero.area_ref; % double, kg, vehicle mass

% Guidance
in.v.gnc.g.p.rate = 1; % nd, guidance rate

in.v.gnc.g.p.bank_mode = uint8(1); %  nd, constant bank
in.v.gnc.g.p.aoa_mode = uint8(1); %  nd, constant angle-of-attack
in.v.gnc.g.p.dm_mode = uint8(6); % manual jettison

% in.v.gnc.g.p_manual.area_ref = [in.v.aero.area_ref; in.v.mp.m_ini/(in.v.aero.cd*beta2)]; % m^2
% in.v.gnc.g.p_manual.mass = in.v.mp.m_ini; % kg

% Control
in.v.gnc.c.p.rate = in.s.traj.rate;

end % manual_earth_in.m
