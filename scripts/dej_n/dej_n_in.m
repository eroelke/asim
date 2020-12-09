% dej_n_in.m 
%   Create input structure for discrete-event jettison aerocapture
%
% Aeroassist Simulation source code (ASIM)
%
% Developed by:
%   Entry systems Design Lab (EsDL)
%   University of Colorado Boulder
%
% Inputs:
%   mode - guidance mode
%   h_atm - atmospheric interface (m)
%   planet - string or index for planet model
%     1: Venus
%     2: Earth
%     3: Mars
%     4: Jupiter
%     5: Saturn
%     6: Titan
%     7: Uranus
%     8: Neptune
%   atm_mode - atmospheric mode
%       1 - exponential atm
%       2 - table look up, no winds
%       3 - table look up, winds
%   
% Outputs:
%   in - struct, multi, input structure with parameters for SEJ guidance
%
% Major Revision History:
%   *Created Jan 2019, E. Roelke

function in = dej_n_in(mode, h_atm,planet,atm_mode)
%#codegen

%% Define input structure
in = define_in;

%% Planetary data
if ischar(planet) == true
    models = {'venus','earth','mars','jupiter','saturn','titan','uranus','neptune'};
    planet = find(strcmpi(planet,models) == 1);
end
if atm_mode > 3 || atm_mode < 1
    fprintf('Warning: Incorrect Atmospheric Mode. Defaulting to Exponential Model\n');
    atm_mode = 1;
end
in.p = get_planetary_data(planet, atm_mode); % Load planet model, GRAM atmosphere + winds


%% Simulation data

% Termination conditions
in.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in.s.term.or.value(1) = h_atm; % double, m, terminate at atmospheric exit
in.s.term.or.direction(1) = int8(1); % int8, nd, positive crossing

in.s.term.or.type(2) = int8(-1); % int8, nd, crossing/equality condition
in.s.term.or.var(2) = uint8(1); % uint8, nd, termination mode (altitude)
in.s.term.or.value(2) = 0.0; % double, m, terminate at atmospheric exit

% Entry State
in.s.traj.gamma_pp = -5*pi/180;
in.s.traj.alt = h_atm; % initial altitude, m
in.s.traj.lat = 0*pi/180; % initial latitude, rad
in.s.traj.lon = 0*pi/180; % initial longitude, rad
in.s.traj.az = 90*pi/180; % initial azimuth, deg->rad
in.s.traj.vel_pp_mag = 10e3;

% Convert scalars to inertial position and velocity vectors
[in.s.traj.r_pci_ini, in.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
    in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
    in.p.omega, in.s.traj.t_ini, 0, in.p.r_e, in.p.r_p);

% Simulation
in.s.traj.rate = 100; % double, Hz, integration rate
in.s.traj.t_ini = 0; % double, s, initial time
in.s.traj.t_max = 500; % double, s, maximum allowable time
in.s.data_rate = 1; % Hz, guidance rate
in.v.gnc.g.p.rate = 1;  % guidance rate
in.v.gnc.c.p.rate = in.s.traj.rate;    % control rate = traj rate
in.v.gnc.n.p.rate = in.s.traj.rate; %nav rate

%% Vehicle data
% Aerodynamics
in.v.aero.mode = uint8(1); % constant aero
in.v.aero.aoa = 0; % double, rad, angle-of-attack
in.v.aero.cl = 0.0; % double, nd, lift coefficient
in.v.aero.cd = 1; % double, nd, drag coefficient
in.v.aero.area_ref = 50; % double, m^2, aerodynamic reference area
in.v.aero.nose_radius = 1; % double, m, nose radius
in.v.mp.m_ini = 1000; % double, kg, vehicle mass

% Guidance
in.v.gnc.g.p.bank_mode = uint8(1); %  nd, constant bank
in.v.gnc.g.p.aoa_mode = uint8(1); %  nd, constant angle-of-attack
in.v.gnc.g.p.prop_mode = uint8(0); %  nd, none
in.v.gnc.g.p.dm_mode = uint8(mode);   % guidance mode

switch (mode)
    case 2 %dcfj
    case 3 %pdg
        in.v.gnc.g.pd.A_sens_atm = 0.5; % m/s2, sensed accel for atmosphere
        in.v.gnc.g.pd.ha_tgt = 400e3;     % m, tgt apoapsis
        in.v.gnc.g.pd.ha_tol = 10e3;   % 50 km tolerancef
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
        
    case 4 %npc
        % DEJ-N settings
        in.v.gnc.g.p_dej_n.npc_mode = uint8(1); % newton method
        in.v.gnc.g.p_dej_n.iter_max = uint8(1); % nd, max NPC internal iterations
        in.v.gnc.g.p_dej_n.A_sens_atm = 0.5; % m/s^2, sensible decel
        in.v.gnc.g.p_dej_n.t_init = 0;  % s, time to start guidance
        in.v.gnc.g.p_dej_n.tgt_ap = 400e3; % m, target apoapsis alt
        in.v.gnc.g.p_dej_n.tol_ap = 10e3; % m, apoapsis alt tolerance
        in.v.gnc.g.p_dej_n.t_max = 700; %s, max predictor integration time step (unused)
        in.v.gnc.g.p_dej_n.n_jett = uint8(1);   % default 1 stage to jettison
        in.v.gnc.g.p_dej_n.tol = 1e-5;  % convergence tolerance
        in.v.gnc.g.p_dej_n.tj_ini = [100 0 0 0 0]'; % s, initial guesses for jettison times
        in.v.gnc.g.p_dej_n.traj_rate = in.s.traj.rate; % trajectory integration rate
        in.v.gnc.g.p_dej_n.area_refs(1:2) = [in.v.aero.area_ref; in.v.aero.area_ref/2]; % m^2
        in.v.gnc.g.p_dej_n.cds(1:2) = in.v.aero.cd; % nd, internal vehicle model cds
        in.v.gnc.g.p_dej_n.cls(1:2) = in.v.aero.cl; % nd, internal vehicle model cls
        in.v.gnc.g.p_dej_n.masses(1:2) = [in.v.mp.m_ini in.v.mp.m_ini/2]; % kg, first two stage masses
    case 6 %manual
        in.v.gnc.g.p_manual.area_refs(1:2) = [in.v.aero.area_ref; in.v.aero.area_ref/2]; % m^2
        in.v.gnc.g.p_manual.masses(1:2) = [in.v.mp.m_ini in.v.mp.m_ini/2]; % kg
        in.v.gnc.g.p_manual.n_jett = uint8(1);
end

% guidance planet model
in.v.gnc.g.p.planet.p_ind = uint8(planet);
in.v.gnc.g.p.planet.r_p = in.p.r_p;
in.v.gnc.g.p.planet.r_e = in.p.r_e;
in.v.gnc.g.p.planet.mu = in.p.mu; % m^3/s^2
in.v.gnc.g.p.planet.j2 = in.p.j2; % nd
in.v.gnc.g.p.planet.npole = [0;0;1]; % nd
in.v.gnc.g.p.planet.omega = in.p.omega; % rad/s
in.v.gnc.g.p.planet.rho0 = in.p.atm.dens_ref;
in.v.gnc.g.p.planet.H = in.p.atm.scale_height;
in.v.gnc.g.p.planet.atm_nom = in.p.atm.table(:,1:7);    %nominal atmosphere
in.v.gnc.g.p.planet.alt_max = h_atm;
in.v.gnc.g.p.planet.alt_min = 0.0;

end % dej_n_earth_in()
