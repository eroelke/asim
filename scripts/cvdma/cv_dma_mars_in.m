% cv_dma_mars_in.m
%   input structure for continuous drag modulation at Mars
% 
% Entry systems Design Lab (EsDL)
% Colorado Center for Astrodynamics Research (CCAR)
% University of Colorado Boulder
% 
% Written By: Evan Roelke, Dec 2018
% 
function in = cv_dma_mars_in()

in = default_in;    % default input structure

% atmospheric data
in.p = get_planetary_data( 3, 3 ); % Load settings for Mars, GRAM atmosphere + winds
atm_data = load('atm_mars_gram2010_1000mc.mat');     % Atmospheric density profiles with dispersion applied
in.p.atm.table = atm_data.nom_table;

%% Simulation data
% Termination conditions
in.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in.s.term.or.value(1) = 125e3; % double, m, terminate at atmospheric exit
in.s.term.or.direction(1) = int8(1); % int8, nd, positive crossing

in.s.term.or.type(2) = int8(-1); % int8, nd, crossing/equality condition
in.s.term.or.var(2) = uint8(1); % uint8, nd, termination mode (altitude)
in.s.term.or.value(2) = 0; % double, m, terminate at atmospheric exit

in.s.traj.rate = 100; % double, Hz, integration rate
in.s.traj.t_ini = 0; % double, s, initial time
in.s.traj.t_max = 500; % double, s, maximum allowable time

% Aerodynamics
in.v.aero.mode = uint8(2); % uint8, nd, aerodynamics mode
in.v.aero.aoa = 0; % double, rad, angle-of-attack
in.v.aero.cl = 0.0; % double, nd, lift coefficient
in.v.aero.cd = 1.53; % double, nd, drag coefficient
in.v.aero.area_ref = 200; % double, m^2, aerodynamic reference area
in.v.aero.nose_radius = 1; % double, m, nose radius
temp = load(['.' filesep 'data' filesep 'aero_sphere_cone_45deg.mat']);
in.v.aero.table = temp.table;

% extract reference trajectory
ref = load('cv_dma/mars/nom_traj.mat');
in.v.aero.cd = ref.cd;
in.v.aero.area_ref = pi*ref.rcs(1)^2;
in.v.mp.m_ini = ref.m(1);

in.s.traj.vel_pp_mag = ref.v_atm*1000;
in.s.traj.gamma_pp = ref.efpa * pi/180;


%% guidance settings
in.v.gnc.g.p.rate = 1; % nd, guidance rate
in.v.gnc.g.p.dm_mode = uint8(5); % continuous drag modulation

% convert reference traj to table form
length = size(ref.out.traj.alt,1);
if length <= 1000
    in.v.gnc.g.p_cvdma.decel(1:length,:) = [ref.out.traj.alt, ref.out.traj.g_loading];
    in.v.gnc.g.p_cvdma.decel(length+1:1000,:) = nan(1000-length,2);
else
% to do
end

% input paremeters
in.v.gnc.g.p_cvdma.mode = uint8(0);
in.v.gnc.g.p_cvdma.iter_max = uint8(1);
in.v.gnc.g.p_cvdma.A_sens_atm = 0.25;
in.v.gnc.g.p_cvdma.ha_tgt = 2000e3;
in.v.gnc.g.p_cvdma.ha_tol = 100e3;
in.v.gnc.g.p_cvdma.cd_ini = in.v.aero.cd;
in.v.gnc.g.p_cvdma.rn = ref.rcs(1)/2;
in.v.gnc.g.p_cvdma.delta_c_ini = 45* pi/180;
in.v.gnc.g.p_cvdma.rc_ini = ref.rcs(1);
in.v.gnc.g.p_cvdma.area_rate = 0.25;
in.v.gnc.g.p_cvdma.rc_bounds = [0.1 2.5];
in.v.gnc.g.p_cvdma.mass = in.v.mp.m_ini;
in.v.gnc.g.p_cvdma.c = ref.rcs(1)/tan(in.v.gnc.g.p_cvdma.delta_c_ini);
in.v.gnc.g.p_cvdma.traj_rate = in.s.traj.rate;


% planet model
in.v.gnc.g.p_cvdma.planet.r_p = in.p.r_p;
in.v.gnc.g.p_cvdma.planet.r_e = in.p.r_e;
in.v.gnc.g.p_cvdma.planet.mu = in.p.mu; % m^3/s^2
in.v.gnc.g.p_cvdma.planet.j2 = in.p.j2; % nd
in.v.gnc.g.p_cvdma.planet.npole = [0;0;1]; % nd
in.v.gnc.g.p_cvdma.planet.omega = in.p.omega; % rad/s
in.v.gnc.g.p_cvdma.planet.atm_table = in.p.atm.table(:,1:7);

    
% overall guidance settings
in.v.gnc.g.p.bank_mode = uint8(1); %  nd, constant bank
in.v.gnc.g.p.aoa_mode = uint8(1); %  nd, constant angle-of-attack
in.v.gnc.g.p.prop_mode = uint8(5); %  nd, SEJ-N guidance
in.v.gnc.g.p.K_bounds = [0.1 2];
in.v.gnc.g.p.K_gain = 0.1;
in.v.gnc.g.p.Kflag = uint8(0);


end