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

function in = dme_mars_guided_in()
%#codegen

%% Define input structure

in = default_in;


%% Planetary data

in.p = get_planetary_data( 3, 3 ); % Load settings for Mars, exponential atmosphere
temp = load([ '.' filesep 'data' filesep 'atm_mars_gram2010_nominal_MSL.mat']);
in.p.atm.table = temp.table;


%% Simulation data

% Termination conditions
in.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in.s.term.or.value(1) = -1500; % double, m, terminate at -1 km altitude
in.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

% in.s.term.or.type(4) = int8(0); % int8, nd, crossing/equality condition
% in.s.term.or.var(4) = uint8(4); % uint8, nd, termination mode (altitude)
% in.s.term.or.value(4) = 500; % double, m, terminate at -2 km MOLA
% in.s.term.or.direction(4) = int8(-1); % int8, nd, negative crossing


% Integration parameters
in.s.traj.rate = 100; % double, Hz, integration rate
in.s.traj.t_ini = 0; % double, s, initial time
in.s.traj.t_max = 1000; % double, s, maximum allowable time
in.s.data_rate = 1;

% Initial state (Eastbound equatorial)
% in.s.traj.alt = 125e3; % initial altitude, m
% in.s.traj.lat = 0; % initial latitude, rad
% in.s.traj.lon = 0; % initial longitude, rad
% in.s.traj.vel_pp_mag = 5500; % initial velocity, m/s
% in.s.traj.gamma_pp = -11*pi/180; % initial flight-path, deg->rad
% in.s.traj.az = 90*pi/180; % initial azimuth, deg->rad
% tgt_lat = 0.0; % rad
% tgt_lon = 930e3/in.p.r_e; % rad

% Initial state and target (Gale Crater, approx.)
in.s.traj.alt = 125e3; % initial altitude, m

% 850 km range
% in.s.traj.lat = -16.125*pi/180; % initial latitude, rad
% in.s.traj.lon = 127.285*pi/180; % initial longitude, rad
% 855 km range
% in.s.traj.lat = -16.188*pi/180; % initial latitude, rad
% in.s.traj.lon = 127.225*pi/180; % initial longitude, rad
% 860 km range
in.s.traj.lat = -16.258*pi/180; % initial latitude, rad
in.s.traj.lon = 127.165*pi/180; % initial longitude, rad
% 865 km range
% in.s.traj.lat = -16.325*pi/180; % initial latitude, rad
% in.s.traj.lon = 127.105*pi/180; % initial longitude, rad

in.s.traj.vel_pp_mag = 5500; % initial velocity, m/s
in.s.traj.gamma_pp = -11*pi/180; % initial flight-path, deg->rad
in.s.traj.az = 45*pi/180; % initial azimuth, deg->rad
tgt_lat = -5.4*pi/180; % rad
tgt_lon = 137.7*pi/180; % rad

% Convert scalars to inertial position and velocity vectors
[in.s.traj.r_pci_ini, in.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
    in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
    in.p.omega, in.s.traj.t_ini, 0, in.p.r_e, in.p.r_p);


%% Vehicle data

% Mass properties
in.v.mp.m_ini = 980; % kg, vehicle mass at entry
in.v.mp.m_jettison = 150; % kg, mass to jettison

% Aerodynamics
in.v.aero.mode = uint8(2); % nd, aerodynamics mode (table look-up)
in.v.aero.area_ref = pi*(4.5/2)^2; % m^2, aerodynamic reference area, assumes 4.7 diameter possible
in.v.aero.nose_radius = 0.66; % m, nose radius: MPF, MER, Phoenix
in.v.aero.cd = 1.68; % nd, from EDL challenges paper, not used
in.v.aero.cl = 0; % nd, not used
% temp = load(['.' filesep 'data' filesep 'aero_sphere_cone_70deg.mat']);
temp = load( ['.' filesep 'data' filesep 'aero_sphere_cone_70deg_aoa_new.mat'] );
cl_zero_aoa = interp1(temp.aoa,temp.cl,0); % interp at zero aoa
cd_zero_aoa = interp1(temp.aoa,temp.cd,0); % interp at zero aoa
mach = linspace(min(temp.mach),max(temp.mach),size(in.v.aero.table,1));
cl_mach = interp1(temp.mach,cl_zero_aoa,mach); % interp at zero aoa
cd_mach = interp1(temp.mach,cd_zero_aoa,mach); % interp at zero aoa
in.v.aero.table(:,1) = mach;
in.v.aero.table(:,2) = cl_mach;
in.v.aero.table(:,3) = cd_mach;

% Decelerator (DGB: MSL aero data, MER sizing)
in.v.decel.mode = uint8(3); % dispersed cd-mach table
in.v.decel.cd = 0.48; % nd, from EDL challenges paper, not used
in.v.decel.d = 14; % m, from EDL challenges paper
in.v.decel.K = 0; % no dispersions (nominal)
temp = load(['.' filesep 'data' filesep 'aero_MSL_DGB.mat']);
in.v.decel.table = temp.table;


%% GNC

% Navigation
in.v.gnc.n.p.mode = uint8(2); % Use Markov process/ECRV error model
in.v.gnc.n.p.rate = 20; % Hz
in.v.gnc.n.p.seed = uint32(1); % nd, Error model seed
in.v.gnc.n.p.tau = 3600; % s, time constant
in.v.gnc.n.p.omega = in.p.omega; % rad/s, planet angular velocity vector
in.v.gnc.n.p.r_e = in.p.r_e; % m, planet equatorial radius
in.v.gnc.n.p.r_p = in.p.r_p; % m, planet polar radius
P_SS = zeros(9);
% P_SS(1:6,1:6) = [ ...
%     2.981933182e-02, 5.369876516e-02, 7.337381389e-02, 3.062525841e-05, 1.173553844e-05, 9.742939383e-06; ...
%     5.369876516e-02, 1.144226269e-01, 1.498759926e-01, 6.301336474e-05, 2.198443978e-05, 1.900429607e-05; ...
%     7.337381389e-02, 1.498759926e-01, 1.983781389e-01, 8.324339215e-05, 2.974134394e-05, 2.542989381e-05; ...
%     3.062525841e-05, 6.301336474e-05, 8.324339215e-05, 3.494835789e-08, 1.243167048e-08, 1.065436605e-08; ...
%     1.173553844e-05, 2.198443978e-05, 2.974134394e-05, 1.243167048e-08, 4.663582047e-09, 3.905447999e-09; ...
%     9.742939383e-06, 1.900429607e-05, 2.542989381e-05, 1.065436605e-08, 3.905447999e-09, 3.309286454e-09 ] * (1000^2);
P_SS(1:6,1:6) = [ ... % MSL estimate from TM/JPL, propagated from EI-9min
   0.552735178453176, -2.302586551541820,  0.912706343673949,  0.000213820256592, -0.001657484328789, -0.000499597478684; ...
  -2.302586551541820,  9.983939761648694, -3.970310667819386, -0.000945977701615,  0.007173593718303,  0.002174779407013; ...
   0.912706343673949, -3.970310667819386,  1.641302368691302,  0.000378135684000, -0.002853771419862, -0.000901624871562; ...
   0.000213820256592, -0.000945977701615,  0.000378135684000,  0.000000095938445, -0.000000679629112, -0.000000207200599; ...
  -0.001657484328789,  0.007173593718303, -0.002853771419862, -0.000000679629112,  0.000005157214989,  0.000001563265207; ...
  -0.000499597478684,  0.002174779407013, -0.000901624871562, -0.000000207200599,  0.000001563265207,  0.000000502604395] * 1.0e+04;
P_SS(7:9,7:9) = eye(3)*(90*in.c.earth_g*1e-6 / 3)^2; % m/s^2, 90 micro-g, 3-sigma
in.v.gnc.n.p.P_SS = P_SS;

% Guidance
in.v.gnc.g.p.rate = 0.5; % Hz
in.v.gnc.g.p.bank_mode = uint8(2); % nd, constant bank rate
in.v.gnc.g.p.const_bank_rate = 4*pi/60; % rad/s, 2 RPM
in.v.gnc.g.p.dma_mode = uint8(7); % nd, SEJ-E guidance
in.v.gnc.g.p.tgt_lat = tgt_lat; % rad, for range calculations
in.v.gnc.g.p.tgt_lon = tgt_lon; % rad, for range calculations
% SEJ-E algorithm
in.v.gnc.g.p_sej_e.mode = uint8(3); % nd, guidance mode
in.v.gnc.g.p_sej_e.iter_max = uint8(5); % nd, maximum iterations
in.v.gnc.g.p_sej_e.A_sens_atm = 0.5; % m/s^2
in.v.gnc.g.p_sej_e.tgt_lat = tgt_lat; % rad
in.v.gnc.g.p_sej_e.tgt_lon = tgt_lon; % rad
% in.v.gnc.g.p_sej_e.tol = 0.5e3; % m, tolerance for corrector
in.v.gnc.g.p_sej_e.tol = 0.1e3; % m, tolerance for corrector
in.v.gnc.g.p_sej_e.cd = 1.7; % nd, average hypersonic CD
in.v.gnc.g.p_sej_e.cl = 0; % nd, average hypersonic CL
in.v.gnc.g.p_sej_e.mass = [980; 830]; % kg, nominal mass before (1) and after (2) skirt jettision
in.v.gnc.g.p_sej_e.flatness = 1 - sqrt(1 - (in.p.r_e^2-in.p.r_p^2)/in.p.r_e^2);
in.v.gnc.g.p_sej_e.v_j_min = 540; % m/s, minimum skirt jettison velocity
in.v.gnc.g.p_sej_e.alt_min = -1500; % m, minimum allowable altitude for predictor
in.v.gnc.g.p_sej_e.alt_max = 130e3; % m, maximum allowable altitude for predictor
in.v.gnc.g.p_sej_e.K_dens_min = 0.2; % nd, minimum allowable density correction factor
in.v.gnc.g.p_sej_e.K_dens_max = 2.0; % nd, maximum allowable density correction factor
% in.v.gnc.g.p_sej_e.K_dens_gain = 0.1; % nd
in.v.gnc.g.p_sej_e.K_dens_gain = 0.05; % nd, density correction factor filter gain
in.v.gnc.g.p_sej_e.t_max = 2000; % s, maximum allowable time for predictor
in.v.gnc.g.p_sej_e.p_r = in.p.r_e; % m, equatorial radius
in.v.gnc.g.p_sej_e.mu = in.p.mu; % m^3/s^2, gravitational parameter
in.v.gnc.g.p_sej_e.j2 = in.p.j2; % nd, j2 perturbation
in.v.gnc.g.p_sej_e.t_jettison_ini = 150; % s, initial jettison time
in.v.gnc.g.p_sej_e.t_inc_max = 1; % s, maximum time increment to adjust jettison time, also nominal predictor step size
in.v.gnc.g.p_sej_e.t_inc_min = 0.05; % s, minimum predictor step size
in.v.gnc.g.p_sej_e.trust_region_ini = 10; % s, maximum allowable change in jettison time when polating
in.v.gnc.g.p_sej_e.area_ref = [pi*(4.5/2)^2; pi*(2.65/2)^2]; % m^2, aero ref area before (1) and after (2) jettison
in.v.gnc.g.p_sej_e.npole = [0;0;1]; % nd, unit vector along north pole of planet
in.v.gnc.g.p_sej_e.omega = in.p.omega; % rad/s, planet rotation rate about npole
in.v.gnc.g.p_sej_e.v_chute = 420; % m/s, nominal parachute deploy velocity
in.v.gnc.g.p_sej_e.chute_area_ref = pi*(in.v.decel.d/2)^2; % m^2, aero ref area for deployed parachute
in.v.gnc.g.p_sej_e.chute_cd = 0.625; % nd, Parachute drag coefficient, based on mach-avg CD, deploy at M=2, touchdown at M=0.3
in.v.gnc.g.p_sej_e.v_chute_box = [420;480]; % m/s, minimum (1) and maximum (2) allowable deploy velocities
% Construct atm table for predictor (based on avg atm density, 100 points only)
% temp = load([ '.' filesep 'data' filesep 'atm_mars_gram2010_nominal_MSL.mat']);
temp = load([ '.' filesep 'data' filesep 'atm_mars_gram2010_average_MSL.mat']);
alt_domain = linspace(min(temp.table(:,1)), 125e3, 100);
for ii = 1:length(alt_domain)
   in.v.gnc.g.p_sej_e.atm_table(ii,1) = alt_domain(ii); % m
   in.v.gnc.g.p_sej_e.atm_table(ii,2) = interp1(temp.table(:,1),temp.table(:,2),alt_domain(ii)); % kg/m^3
end
clear temp;

% Control
in.v.gnc.c.p.rate = 20; % Hz
in.v.gnc.c.p.decel_mode = uint8(1);

end % dme_mars_in()