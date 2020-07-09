


function in = pred_guid_venus_in()
%% Define Input Structure
in = define_in;
%% Planetary Data
in.p = get_planetary_data(1,3); % venus, with winds
temp = load('.\data\atm_data\atm_venus_mc.mat'); % nominal VenusGRAM table from JPL
in.p.atm.table = temp.nom_table;

%% Termination Conditions
in.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in.s.term.or.value(1) = 150e3; % double, m, terminate at atmospheric exit
in.s.term.or.direction(1) = int8(1); % int8, nd, positive crossing

in.s.term.or.type(2) = int8(-1); % int8, nd, crossing/equality condition
in.s.term.or.var(2) = uint8(1); % uint8, nd, termination mode (altitude)
in.s.term.or.value(2) = 0.0; % double, m, terminate at atmospheric exit

%% Simulation Parameters
in.s.traj.rate = 100; % double, Hz, integration rate
in.s.traj.t_ini = 0; % double, s, initial time
in.s.traj.t_max = 1000; % double, s, maximum allowable time
in.s.data_rate = in.s.traj.rate;


%% Entry State
in.s.traj.gamma_pp = -5.5*pi/180;   % entry flight path angle, rad
in.s.traj.alt = 150e3; % initial altitude, m
% in.s.traj.r_pci_ini = [(in.p.r_e + in.s.traj.alt); 0; 0]; % double(3), m, initial PCI position vector
in.s.traj.lat = 0*pi/180; % initial latitude, rad
in.s.traj.lon = 0*pi/180; % initial longitude, rad
in.s.traj.az = 90*pi/180; % initial azimuth, deg->rad
in.s.traj.vel_pp_mag = 10.5e3;  %m/s, entry velocity

%% Convert scalars to inertial position and velocity vectors
[in.s.traj.r_pci_ini, in.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
    in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
    in.p.omega, in.s.traj.t_ini, 0, in.p.r_e, in.p.r_p);

%% Vehicle Data
in.v.mp.m_ini = 400;    %kg, initial mass
in.v.aero.area_ref = 200; % double, m^2, aerodynamic reference area
in.v.aero.nose_radius = 1; % double, m, nose radius
% Default Aerodynamics
in.v.aero.mode = uint8(1); % uint8, nd, aerodynamics mode (mach-dependent table)
in.v.aero.aoa = 0; % double, rad, angle-of-attack (ballistic)
in.v.aero.cl = 0.0; % double, nd, lift coefficient
in.v.aero.cd = 1.53; % double, nd, drag coefficient
temp = load('./data/aero_data/aero_sphere_cone_45deg.mat');
in.v.aero.table = temp.table;

%% Trajectory Rates
in.v.gnc.g.p.rate = in.s.traj.rate; %guidance rate (Hz)
in.v.gnc.n.p.rate = in.s.traj.rate;
in.v.gnc.c.p.rate = in.s.traj.rate;

%% Guidance Basics
in.v.gnc.g.p.bank_mode = uint8(1);  % constant bank
in.v.gnc.g.p.aoa_mode = uint8(1);   % constant aoa
% in.v.gnc.g.p.prop_mode = uint8(12); % predictive guidance algorithm
in.v.gnc.g.p.dm_mode = uint8(3);

% predictor algorithm data struct
in.v.gnc.g.pd.A_sens_atm = 0.5; % m/s2, sensed accel for atmosphere
in.v.gnc.g.pd.ha_tgt = 10000e3;     % m, tgt apoapsis
in.v.gnc.g.pd.ha_tol = 50e3;   % 50 km tolerancef
in.v.gnc.g.pd.cds = in.v.aero.cd.*ones(6,1);
in.v.gnc.g.pd.cls = in.v.aero.cl.*ones(6,1);
in.v.gnc.g.pd.n_jett = uint8(1);     % single jettison event
in.v.gnc.g.pd.masses(1:2) = [in.v.mp.m_ini 200]';    % spacecraft masses, kg
in.v.gnc.g.pd.mcflag = uint8(0);
in.v.gnc.g.pd.t_max = in.s.traj.t_max;   % max predictor integration time, s
in.v.gnc.g.pd.tj0 = [90 0 0 0 0]';  % initial guess for jettison times, s
in.v.gnc.g.pd.rate = in.s.traj.rate;
in.v.gnc.g.pd.trigger = uint8(0);
in.v.gnc.g.pd.area_refs(1:2) = [in.v.aero.area_ref in.v.aero.area_ref/2]';

in.v.gnc.g.p.planet.alt_max = in.s.traj.alt;
in.v.gnc.g.p.planet.alt_min = 0;
in.v.gnc.g.p.planet.mu = in.p.mu;
in.v.gnc.g.p.planet.npole = [0;0;1];
in.v.gnc.g.p.planet.omega = in.p.omega;
in.v.gnc.g.p.planet.r_p = in.p.r_p;
in.v.gnc.g.p.planet.r_e = in.p.r_e;
in.v.gnc.g.p.planet.j2 = in.p.j2;
in.v.gnc.g.p.planet.atm_nom = in.p.atm.table;

end