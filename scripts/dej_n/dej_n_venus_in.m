% dej_n_earth_in.m 
%   Create input structure for single-event jettison guided aerocapture at
%   Venus
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
%   *Created Sep 2018, E. Roelke

function in = dej_n_venus_in()
%#codegen

%% Define input structure
in = define_in;

%% Planetary data
in.p = get_planetary_data( 1, 3 ); % Load settings for Venus, GRAM atmosphere + winds


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
in.s.traj.t_max = 700; % double, s, maximum allowable time
in.s.traj.r_pci_ini = [(in.p.r_e + alt0); 0; 0]; % double(3), m, initial PCI position vector, 125 km altitude
in.s.traj.v_pci_ini = v0*[sind(fpa0); cosd(fpa0); 0]; % double(3), m/s, initial PCI inertial velocity vector

in.s.data_rate = 1; % Hz, guidance rate

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
temp = load('./data/aero_data/aero_sphere_cone_45deg.mat');
in.v.aero.table = temp.table;

in.v.mp.m_ini = beta1*in.v.aero.cd*in.v.aero.area_ref; % double, kg, vehicle mass

% Guidance
in.v.gnc.g.p.rate = 1; % nd, guidance rate

in.v.gnc.g.p.bank_mode = uint8(1); %  nd, constant bank
in.v.gnc.g.p.aoa_mode = uint8(1); %  nd, constant angle-of-attack
in.v.gnc.g.p.dm_mode = uint8(4); %  nd, DEJ-N guidance

% DEJ-N settings
in.v.gnc.g.p_dej_n.npc_mode = uint8(0);
in.v.gnc.g.p_dej_n.iter_max = uint8(5); % nd
in.v.gnc.g.p_dej_n.A_sens_atm = 0.25; % m/s^2
in.v.gnc.g.p_dej_n.tgt_ap = 400e3; % m
in.v.gnc.g.p_dej_n.tol_ap = 10e3; % m
in.v.gnc.g.p_dej_n.planet.alt_min = 0; % m
in.v.gnc.g.p_dej_n.planet.alt_max = 150e3; % m
% in.v.gnc.g.p_dej_n.K_gain = 0.1;
% in.v.gnc.g.p_dej_n.K_bounds = [0.1, 2.0]; % nd
in.v.gnc.g.p_dej_n.n_jett = uint8(1);   % default 1 stage to jettison
in.v.gnc.g.p_dej_n.tol = 1e-5;  % convergence tolerance
in.v.gnc.g.p_dej_n.tj_ini = [100 0 0 0 0]'; % s, initial guesses for jettison times
in.v.gnc.g.p_dej_n.traj_rate = in.s.traj.rate; % trajectory integration rate
in.v.gnc.g.p_dej_n.area_refs(1:2) = [in.v.aero.area_ref; in.v.mp.m_ini/(in.v.aero.cd*beta2)]; % m^2
in.v.gnc.g.p_dej_n.cds(1) = in.v.aero.cd; % nd
in.v.gnc.g.p_dej_n.cls(1) = in.v.aero.cl; % nd
in.v.gnc.g.p_dej_n.masses(1:2) = [in.v.mp.m_ini in.v.mp.m_ini/2]; % kg

% planet model
in.v.gnc.g.p_dej_n.planet.r_p = in.p.r_p;
in.v.gnc.g.p_dej_n.planet.r_e = in.p.r_e;
in.v.gnc.g.p_dej_n.planet.mu = in.p.mu; % m^3/s^2
in.v.gnc.g.p_dej_n.planet.j2 = in.p.j2; % nd
in.v.gnc.g.p_dej_n.planet.npole = [0;0;1]; % nd
in.v.gnc.g.p_dej_n.planet.omega = in.p.omega; % rad/s
in.v.gnc.g.p_dej_n.planet.atm_table = in.p.atm.table(:,1:7);

% Control
in.v.gnc.c.p.rate = in.s.traj.rate;

end % dej_n_earth_in()
