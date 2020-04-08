% sej_a_earth_in.m 
%   Create input structure for single-event jettison guided aerocapture at
%   Earth
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
%   in - struct, multi, input structure with parameters for SEJ guidance
%
% Major Revision History:
%   *Created NOV 2011, Z.R. Putnam
%   *1 AUG 2012, Z.R.Putnam, updated for inclusion in repository

function in = sej_a_earth_in()
%#codegen

%% Define input structure
in = define_in;


%% Planetary data

in.p = get_planetary_data( 2, 3 ); % Load settings for Earth, table atm model

temp = load('C:\Users\Evan Roelke\Documents\research\asim\data\atm_earth_gram2007_nominal.mat');
in.p.atm.table = temp.table;

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
temp = load(['.' filesep 'data' filesep 'aero_sphere_cone_60deg.mat']);
in.v.aero.table = temp.table;

in.v.mp.m_ini = beta1*in.v.aero.cd*in.v.aero.area_ref; % double, kg, vehicle mass

% Guidance
in.v.gnc.g.p.rate = 1; % nd, guidance rate

in.v.gnc.g.p.bank_mode = uint8(1); %  nd, constant bank
in.v.gnc.g.p.aoa_mode = uint8(1); %  nd, constant angle-of-attack
in.v.gnc.g.p.prop_mode = uint8(0); %  nd, SEJ-A guidance
in.v.gnc.g.p.dm_mode = uint8(1);

% in.v.gnc.g.p_sej_a.mode = uint8(1);
% in.v.gnc.g.p_sej_a.iter_max = uint8(5); % nd
% in.v.gnc.g.p_sej_a.A_sens_atm = 0.25; % m/s^2
% in.v.gnc.g.p_sej_a.tgt_ap = 400e3; % m
% in.v.gnc.g.p_sej_a.tol_ap = 10e3; % m
% in.v.gnc.g.p_sej_a.alt_min = 0; % m
% in.v.gnc.g.p_sej_a.alt_max = 125e3; % m
% in.v.gnc.g.p_sej_a.K_dens_min = 0.1; % nd
% in.v.gnc.g.p_sej_a.K_dens_max = 2.0; % nd
% in.v.gnc.g.p_sej_a.K_dens_gain = 0.1; % nd
% in.v.gnc.g.p_sej_a.t_max = 1500; % s
% in.v.gnc.g.p_sej_a.tj_ini = [100; 1000]; % s
% in.v.gnc.g.p_sej_a.t_inc_max = 2.0; % s
% in.v.gnc.g.p_sej_a.t_inc_min = 0.01; % s
% in.v.gnc.g.p_sej_a.area_ref = [in.v.aero.area_ref; in.v.mp.m_ini/(in.v.aero.cd*beta2)]; % m^2
% in.v.gnc.g.p_sej_a.cd = in.v.aero.cd; % nd
% in.v.gnc.g.p_sej_a.cl = in.v.aero.cl; % nd
% in.v.gnc.g.p_sej_a.mass = in.v.mp.m_ini; % kg
% in.v.gnc.g.p_sej_a.p_r = in.p.r_e; % m
% in.v.gnc.g.p_sej_a.mu = in.p.mu; % m^3/s^2
% in.v.gnc.g.p_sej_a.j2 = in.p.j2; % nd
% in.v.gnc.g.p_sej_a.npole = [0;0;1]; % nd
% in.v.gnc.g.p_sej_a.omega = in.p.omega; % rad/s
% in.v.gnc.g.p_sej_a.atm_table = in.p.atm.table(:,1:7);
% 
% in.v.gnc.g.p_sej_a.trust_region_ini = 10;

% Control
in.v.gnc.c.p.rate = in.s.traj.rate;

end % sej_a_earth_in()
