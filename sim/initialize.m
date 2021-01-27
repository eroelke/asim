% initialize.m 
%   This function is used to initialize all the matrices used in the Aeroassist 
%   simulator. This initialization must be performed if compiling the code to C.
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   in - data structure, misc, input data structure
%   
% Outputs:
%   dat - data strucutre, misc, initialized trajectory data structure (NaNed to
%   length of the time vector) 
%
% Major Revision History:
%   *Created 2008, M. Grant
%   *7 NOV 2011, Z.R. Putnam, added seperate data rate independent of sim
%       integration rate, removed t_vec from dat structure
%   *23 NOV 2011, Z.R. Putnam, Removed unused structures, added termination
%       structure
%   *12 DEC 2011, Z.R. Putnam, added MSL guidance states to output
%   *8 FEB 2012, Z.R. Putnam, added angle-of-attack and command to output
%   *Jul 2018, pre-defined data sizes, E. Roelke
%   *Feb 2019, changed guidance states to reduce overall data output, E. Roelke

function [dat, veh] = initialize(in) %#codegen

%% Vehicle data structure
veh = define_veh;
% Parameters
% veh.p.m_jettison = in.v.mp.m_jettison;
% Non-integrated states
veh.s.not_jettisoned = true; % no jettison event has occured
veh.s.compute_trim = uint8(1); % enable static trim calculation
veh.s.bank = in.s.traj.bank_ini; % rad
veh.s.aoa = in.s.traj.aoa_ini; % rad
veh.s.ssa = in.s.traj.ssa_ini; % rad
veh.s.area_ref = in.v.aero.area_ref; % m^2
veh.s.mass = in.v.mp.m_ini; %initial mass


%% Output Data storage
% Determine storage array size
data_t_step = 1/in.s.data_rate;
t_vec = (in.s.traj.t_ini : data_t_step : in.s.traj.t_max)';
if (in.v.gnc.g.p.dm_mode > 0)
    t_length = length(t_vec) + double(in.v.gnc.g.p_dej_n.n_jett); % +1 for jettison and terminal point
else
    t_length = length(t_vec) + 1;   % +1 for terminal point
end
ldat1 = nan(t_length,1); %1D arrays
% ldat2 = nan(t_length,2);  %2D arrays. unused
ldat3 = nan(t_length,3); %3D arrays

% Trajectory
dat.traj.pos_ii = ldat3; % Inertial position
% dat.traj.pos_ii_mag = ldat1; % Inertial position magnitude
dat.traj.vel_ii = ldat3; % Inertial velocity
dat.traj.vel_ii_mag = ldat1; % Inertial velocity magnitude
dat.traj.alt = ldat1; % Altitude
dat.traj.rho = ldat1; % Density
dat.traj.mass = ldat1; % Mass
dat.traj.force_ii = ldat3; % Inertial external force vector
% dat.traj.thrust_ii = ldat3; % Total inertial thrust vector
% dat.traj.thrust_pp = ldat3; % Total planet relative thrust vector
% dat.traj.thrust_pp_mag = ldat1; % Total planet relative thrust vector magnitude
dat.traj.pos_pp = ldat3; % Planet relative position
% dat.traj.pos_pp_mag = ldat1; % Planet relative position magnitude
dat.traj.V_pf_pci = ldat3; % Planet-fixed velocity
dat.traj.vel_pp = ldat3; % Planet relative velocity
dat.traj.vel_pp_mag = ldat1; % Planet relative velocity magnitude
% dat.traj.gamma_ii = ldat1; % Inertial flight path angle
dat.traj.gamma_pp = ldat1; % Planet relative flight path angle
dat.traj.g_loading = ldat1; % G-loading
% dat.traj.mach = ldat1; % Mach, nd
dat.traj.aoa = ldat1; % angle-of-attack, rad
dat.traj.bank = ldat1; % bank, rad
% dat.traj.bank_rate = ldat1; % bank rate, rad/s
% dat.traj.bank_accel = ldat1; % bank acceleration, rad/s/s
% dat.traj.ssa = ldat1; % sideslip angle, rad
% dat.traj.cmd_aoa = ldat1; % commanded angle-of-attack, rad
% dat.traj.cmd_bank = ldat1; % commanded bank, rad
dat.traj.cl = ldat1; % lift coefficient
dat.traj.cd = ldat1; % drag coefficient
% dat.traj.cd_decel = ldat1; % drag coefficient
% dat.traj.LoD = ldat1; % L/D
dat.traj.gravity_ii = ldat3;
% dat.traj.drag_ii = ldat3;
% dat.traj.lift_ii = ldat3;
dat.traj.drag_ii_mag = ldat1;
dat.traj.lift_ii_mag = ldat1;
% dat.traj.decel_drag_ii = ldat3;
% dat.traj.q = ldat1;
% dat.traj.azi_pp = ldat1;
% dat.traj.downrange = ldat1;
% dat.traj.crossrange = ldat1;
% dat.traj.heat_rate = ldat1;
% dat.traj.lat = ldat1;
% dat.traj.lon = ldat1;
% dat.traj.cg = ldat3;
dat.traj.time = ldat1; % Time
% dat.traj.prop_mass = ldat1;
% dat.traj.h_ap = ldat1;   % way to visualize aerocapture energy loss

% general guidance params
dat.g.K_dens = ldat1;   % density corrector
dat.g.K_true = ldat1;   % true density variation
dat.g.rho_est = ldat1;  % true density estimate
dat.g.atm_hist = nan(1000,4); % internal atmosphere model
dat.g.rss_nom = nan(1,1);  % nominal atm rss error (static)
dat.g.rss_K = ldat1;  % density factor rss err
dat.g.rss_ens = ldat1;    % ensemble rss err
dat.g.ind_curr = ldat1;     % ensemble current mc index

% Vehicle
dat.veh.area_ref = ldat1; % m^2


% Navigation
dat.nav.r_pci = ldat3;
dat.nav.v_inrtl_pci = ldat3;
dat.nav.a_sens_pci = ldat3;
dat.nav.rva_error = nan(t_length,9);
% dat.nav.ercv = nan(t_length, 9);
dat.nav.x_ercv = nan(t_length, 9);
% dat.nav.r_pcpf = ldat3;
% dat.nav.v_pf_pci = ldat3;
% dat.nav.v_pf_pcpf = ldat3;
% dat.nav.a_sens_pcpf = ldat3;
% dat.nav.lat_gd = ldat1;
% dat.nav.lon = ldat1;
% dat.nav.alt = ldat1;
% dat.nav.dt = ldat1;
% dat.nav.t = ldat1;
% dat.nav.t_ini = ldat1;

% Control
% dat.ctrl.deploy_decel = ldat1;
% dat.ctrl.jettison = ldat1;

% guidance data allocation
%     broken into switch statements to reduce amount of data storage
%     E. Roelke, Feb 2019
%% bank angle guidance init
switch in.v.gnc.g.p.bank_mode
    case 4 %msl guidance
        dat.g.msl.tgt_overflown = ldat1;
        dat.g.msl.phase = ldat1;
        dat.g.msl.next_step = ldat1;
        dat.g.msl.cmd_bank_sign = ldat1;
        dat.g.msl.alt_rate = ldat1;
        dat.g.msl.ang_to_tgt = ldat1;
        dat.g.msl.cmd_bank = ldat1;
        dat.g.msl.cmd_lodv = ldat1;
        dat.g.msl.D_mag = ldat1;
        dat.g.msl.lat_ang = ldat1;
        dat.g.msl.lat_db = ldat1;
        dat.g.msl.lod_work = ldat1;
        dat.g.msl.lodv_lim_work = ldat1;
        dat.g.msl.range_final = ldat1;
        dat.g.msl.range_to_tgt = ldat1;
        dat.g.msl.ref_alt_rate = ldat1;
        dat.g.msl.ref_D_mag = ldat1;
        dat.g.msl.ref_f1 = ldat1;
        dat.g.msl.ref_f2 = ldat1;
        dat.g.msl.ref_f3 = ldat1;
        dat.g.msl.ref_range = ldat1;
        dat.g.msl.tgt_transfer_ang = ldat1;
        dat.g.msl.time = ldat1;
        dat.g.msl.time_init = ldat1;
        dat.g.msl.V_mag = ldat1;
        dat.g.msl.V_norm_sq = ldat1;
        dat.g.msl.U_E_pole = ldat3;
        %             dat.g.msl.U_norm_traj(i,:) = ldat3;
        dat.g.msl.U_tgt = ldat3;
        dat.g.msl.U_tgt_init = ldat3;
        dat.g.msl.U_tgt_east = ldat3;
        dat.g.msl.U_tgt_eq = ldat3;
        dat.g.msl.V_work = ldat3;
        
    otherwise
        dat.g.msl = msl_empty();
end %switch bank_mode

%% prop mode guidance init
switch in.v.gnc.g.p.prop_mode
    
    case 8 %adm guidance
        dat.g.adm.jettison = ldat1;
        dat.g.adm.phase = ldat1;
        dat.g.adm.next_step = ldat1;
        dat.g.adm.iter = ldat1;
        dat.g.adm.A_mag = ldat1;
%         dat.g.adm.K_dens = ldat1;
        dat.g.adm.dens_est = ldat1;
        dat.g.adm.V_pf_mag = ldat1;
        dat.g.adm.t_inc = ldat1;
        dat.g.adm.trust_region = ldat1;
        dat.g.adm.ngood = ldat1;
        dat.g.adm.t_j_curr = ldat1;
        dat.g.adm.range_to_tgt = ldat1;
        dat.g.adm.flight_range = ldat1;
        dat.g.adm.delta_range = ldat1;
        dat.g.adm.V_chute = ldat1;
        dat.g.adm.deploy_decel = ldat1;
        dat.g.adm.gs1 = ldat1;
        dat.g.adm.gs2 = ldat1;
        
    otherwise
        dat.g.adm = adm_empty();
end %switch prop mode

%% drag modulation guidance init
switch in.v.gnc.g.p.dm_mode
    case 2 %dcfj
        dat.g.dcfj.stage = ldat1;
        dat.g.dcfj.phase = ldat1;
        dat.g.dcfj.next_step = ldat1;
        dat.g.dcfj.A_mag = ldat1;
        dat.g.dcfj.tj_curr = ldat1;
        dat.g.dcfj.g1_time = ldat1;
        dat.g.dcfj.g2_time = ldat1;
        dat.g.dcfj.g2 = ldat1;
        dat.g.dcfj.togo = ldat1;
        dat.g.dcfj.g_ratio = ldat1;
        
        dat.g.pd = pd_empty();
        dat.g.dej_n = dej_n_empty();
        dat.g.cvdma = cvdma_empty();
        dat.g.manual = manual_empty();
        dat.g.sej_e = sej_e_empty();
        
    case 3 %predictive (pd)
        dat.g.dcfj = dcfj_empty();
        
        dat.g.pd.jettison = nan(t_length,5);
        dat.g.pd.A_mag = ldat1;
%         dat.g.pd.K_dens = ldat1;
        dat.g.pd.ha_j = ldat1;
        dat.g.pd.delta_ap = ldat1;
        dat.g.pd.stage = ldat1;
        dat.g.pd.trig_val = ldat1;
        dat.g.pd.dr_ap = ldat1;
        dat.g.pd.dr_ap_true = ldat1;
        
        dat.g.dej_n = dej_n_empty();
        dat.g.cvdma = cvdma_empty();
        dat.g.manual = manual_empty();
        dat.g.sej_e = sej_e_empty();
        
    case 4  % DEJ-n guidance
        dat.g.dcfj = dcfj_empty();
        dat.g.pd = pd_empty();
        
        dat.g.dej_n.stage = ldat1;
        dat.g.dej_n.iter = ldat1;
        dat.g.dej_n.A_mag = ldat1;
%         dat.g.dej_n.K_dens = ldat1;
        dat.g.dej_n.V_pf_mag = ldat1;
        dat.g.dej_n.t_inc = ldat1;
        dat.g.dej_n.tj = ldat1;
        dat.g.dej_n.r_ap = ldat1;
        dat.g.dej_n.dr_ap = ldat1;
        dat.g.dej_n.dr_ap_true = ldat1;
        dat.g.dej_n.step_size = ldat1;
        dat.g.dej_n.dr_f = nan(5,1);
        dat.g.dej_n.ncalls = ldat1;
        
        %dat.g.dej_n.atm_inds = zeros(1000,1000);
        dat.g.dej_n.dtj = ldat1;
        dat.g.dej_n.dtj_lim = ldat1;

        dat.g.cvdma = cvdma_empty();
        dat.g.manual = manual_empty();
        dat.g.sej_e = sej_e_empty();
        
    case 5  %cv dma guidance
        dat.g.dcfj = dcfj_empty();
        dat.g.pd = pd_empty();
        dat.g.dej_n = dej_n_empty();
        
        dat.g.cvdma.cd = ldat1;
        dat.g.cvdma.Aref = ldat1;
%         dat.g.cvdma.K_dens = ldat1;
        dat.g.cvdma.rc = ldat1;
        dat.g.cvdma.delta_c = ldat1;
        dat.g.cvdma.ra = ldat1;
        dat.g.cvdma.delta_ra = ldat1;
        dat.g.cvdma.beta = ldat1;
        
        dat.g.manual = manual_empty();
        dat.g.sej_e = sej_e_empty();
        
    case 6 % manual
        dat.g.dcfj = dcfj_empty();
        dat.g.pd = pd_empty();
        dat.g.dej_n = dej_n_empty();
        dat.g.cvdma = cvdma_empty();
        
        dat.g.manual.stage = ldat1;
        dat.g.manual.mass = ldat1;
        dat.g.manual.area = ldat1;
        
        dat.g.sej_e = sej_e_empty();
        
    case 7 %sej-e
        dat.g.dcfj = dcfj_empty();
        dat.g.pd = pd_empty();
        dat.g.dej_n = dej_n_empty();
        dat.g.cvdma = cvdma_empty();
        dat.g.manual = manual_empty();
        
        dat.g.sej_e.jettison = ldat1;
        dat.g.sej_e.phase = ldat1;
        dat.g.sej_e.next_step = ldat1;
        dat.g.sej_e.iter = ldat1;
        dat.g.sej_e.A_mag = ldat1;
%         dat.g.sej_e.K_dens = ldat1;
        dat.g.sej_e.dens_est = ldat1;
        dat.g.sej_e.V_pf_mag = ldat1;
        dat.g.sej_e.t_inc = ldat1;
        dat.g.sej_e.trust_region = ldat1;
        dat.g.sej_e.ngood = ldat1;
        dat.g.sej_e.t_j_curr = ldat1;
        dat.g.sej_e.range_to_tgt = ldat1;
        dat.g.sej_e.flight_range = ldat1;
        dat.g.sej_e.delta_range = ldat1;
        dat.g.sej_e.v_chute = ldat1;
        dat.g.sej_e.deploy_decel = ldat1;
        
    otherwise %no guidance
        dat.g.dcfj = dcfj_empty();
        dat.g.pd = pd_empty();
        dat.g.dej_n = dej_n_empty();
        dat.g.cvdma = cvdma_empty();
        dat.g.manual = manual_empty();
        dat.g.sej_e = sej_e_empty();
end %siwtch dm_mode


% Termination
dat.term.cond = uint8(zeros(10,1));

% SNC MC data
% dat.snc.alt = ldat1; % altitude, ft
% dat.snc.fpa = ldat1; % flight-path angle, deg
% dat.snc.dynp = ldat1; % dynamic pressure, psf
% dat.snc.mach = ldat1; % mach number, nd
% dat.snc.v_wind_mag = ldat1; % air-relative velocity magnitude
% dat.snc.v_inrtl = ldat3; % inertial velocity vector
% dat.snc.azi = ldat1; % azimuth, deg east of north
% dat.snc.lon = ldat1; % longitude, deg
% dat.snc.lat = ldat1; % geodetic latitude, deg
% dat.snc.weight = ldat1; % weight vector magnitude, lbf
% dat.snc.alpha = ldat1; % angle-of-attack, deg
% dat.snc.beta = ldat1; % sideslip angle, deg
% dat.snc.bank = ldat1; % bank angle, deg
% dat.snc.range = ldat1; % downrange to tgt, from SNC definition, nmi
% dat.snc.crossrange = ldat1; % crossrange to tgt, from SNC definition, nmi
% dat.snc.cg = ldat3; % c.g. location, in
% dat.snc.bank_cmd = ldat1; % bank command, deg
% dat.snc.alt_rate = ldat1; % altitude rate, ft/s
% dat.snc.heat_rate = ldat1; % 1-ft ref. sphere stagnation point convective heating, BTU/ft^2/s
% dat.snc.time = ldat1; % reference time from EI, s
% dat.snc.ground_speed = ldat1; % ground speed, fps
% dat.snc.pres = ldat1; % ambient pressure, psf
% dat.snc.temp = ldat1; % ambient temperature, R
% dat.snc.dens = ldat1; % ambient density, lbm/in3
% dat.snc.wind = ldat3; % winds, ft/s

end % initialize()


function msl = msl_empty()
msl.tgt_overflown = nan;
msl.phase = nan;
msl.next_step = nan;
msl.cmd_bank_sign = nan;
msl.alt_rate = nan;
msl.ang_to_tgt = nan;
msl.cmd_bank = nan;
msl.cmd_lodv = nan;
msl.D_mag = nan;
msl.lat_ang = nan;
msl.lat_db = nan;
msl.lod_work = nan;
msl.lodv_lim_work = nan;
msl.range_final = nan;
msl.range_to_tgt = nan;
msl.ref_alt_rate = nan;
msl.ref_D_mag = nan;
msl.ref_f1 = nan;
msl.ref_f2 = nan;
msl.ref_f3 = nan;
msl.ref_range = nan;
msl.tgt_transfer_ang = nan;
msl.time = nan;
msl.time_init = nan;
msl.V_mag = nan;
msl.V_norm_sq = nan;
msl.U_E_pole = nan;
%             msl.U_norm_traj(i,:) = ldat3;
msl.U_tgt = nan;
msl.U_tgt_init = nan;
msl.U_tgt_east = nan;
msl.U_tgt_eq = nan;
msl.V_work = nan;    
end

function adm = adm_empty()
adm.jettison = nan;
adm.phase = nan;
adm.next_step = nan;
adm.iter = nan;
adm.A_mag = nan;
% adm.K_dens = nan;
adm.dens_est = nan;
adm.V_pf_mag = nan;
adm.t_inc = nan;
adm.trust_region = nan;
adm.ngood = nan;
adm.t_j_curr = nan;
adm.range_to_tgt = nan;
adm.flight_range = nan;
adm.delta_range = nan;
adm.V_chute = nan;
adm.deploy_decel = nan;
adm.gs1 = nan;
adm.gs2 = nan;
end

function dcfj = dcfj_empty()
dcfj.stage = nan;
dcfj.phase = nan;
dcfj.next_step = nan;
dcfj.A_mag = nan;
dcfj.tj_curr = nan;
dcfj.g1_time = nan;
dcfj.g2_time = nan;
dcfj.g2 = nan;
dcfj.togo = nan;
dcfj.g_ratio = nan;
end

function pd = pd_empty()
pd.jettison = nan;
pd.A_mag = nan;
% pd.K_dens = nan;
pd.ha_j = nan;
pd.delta_ap = nan;
pd.stage = nan;
pd.trig_val = nan;
pd.dr_ap = nan;
pd.dr_ap_true = nan;
end

function dej_n = dej_n_empty()
dej_n.stage = nan;
dej_n.iter = nan;
dej_n.A_mag = nan;
% dej_n.K_dens = nan;
dej_n.V_pf_mag = nan;
dej_n.t_inc = nan;
dej_n.tj = nan;
dej_n.r_ap = nan;
dej_n.dr_ap = nan;
dej_n.dr_ap_true = nan;
dej_n.step_size = nan;
dej_n.dr_f = nan;
dej_n.ncalls = nan;
dej_n.dtj = nan;
dej_n.dtj_lim = nan;
end

function cvdma = cvdma_empty()
cvdma.cd = nan;
cvdma.Aref = nan;
% cvdma.K_dens = nan;
cvdma.rc = nan;
cvdma.delta_c = nan;
cvdma.ra = nan;
cvdma.delta_ra = nan;
cvdma.beta = nan;
end

function manual = manual_empty()
manual.stage = nan;
manual.mass = nan;
manual.area = nan;
end

function sej_e = sej_e_empty()
sej_e.jettison = nan;
sej_e.phase = nan;
sej_e.next_step = nan;
sej_e.iter = nan;
sej_e.A_mag = nan;
% sej_e.K_dens = nan;
sej_e.dens_est = nan;
sej_e.V_pf_mag = nan;
sej_e.t_inc = nan;
sej_e.trust_region = nan;
sej_e.ngood = nan;
sej_e.t_j_curr = nan;
sej_e.range_to_tgt = nan;
sej_e.flight_range = nan;
sej_e.delta_range = nan;
sej_e.v_chute = nan;
sej_e.deploy_decel = nan;
end

