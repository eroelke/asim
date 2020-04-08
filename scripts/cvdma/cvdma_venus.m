%% cv_dma_venus
%   -Venus drag modulation aerocapture trajectory.
%   -continous drag modulation control
%
% Developed by:
%   Entry systems Design Lab (EsDL)
%       Colorado Center for Astrodynamics Research
%       CU Boulder
% 
% Inputs:
%   x0: entry state struct
%       efpa: entry flight-path angle (deg)
%       v0: planet-relative entry speed (km/s)
%       h0: entry altitude (km)
%       lat0: entry latitude (deg)
%       lon0: entry longitude (deg)
%       az0: entry azimuth (deg)
%   gnc: guidance struct
%       ha_tgt: target apoapsis altitude (km)
%       ha_tol: apoapsis altitude tolerance (km)
%       area_rate: max/min area change rate (m per guidance step)
%       rc_bounds: [min max] backshell radii values
%       guid_mode:
%           0: reference trajectory mode (incomplete)
%           1: numerical predictor-corrector (NPC)
%       guid_rate: guidance call rate (Hz)
%   aero: geometry and aerodynamics struct
%       rc0: entry backshell radius (m)
%       m: mass (kg) - assumed constant throughout traj
%       cd0: initial drag coefficient
%       delta_c0: initial half-cone angle (deg)
%       rn: nose radius (m) - assumed constant
%   sim: simulation inputs data struct
%       traj_rate: trajectory integration rate (Hz)
%       data_rate: data storage rate (Hz)
%   mc: monte carlo input struct
%       mc.flag:    flag for monte carlo run (and density corrector)
%       mc.N:       number of monte carlo runs
%       mc.sigs:    [m (kg),cl,cd,aoa,v-atm (m/s),h_atm (m),fpa (rad),lat,lon,az]
% 
% Outputs:
%   out:        simulation output data structure
%
% Major Revision History:
%   Feb 2019 - Created, Evan Roelke
% 
function out = cvdma_venus(x0,aero,gnc,sim,mc)
format longg

%% Geometry and Aerodynamics
area_ini = pi*aero.rc0^2;
% rn = rc_ini/2;     % Nose radius (m)

% chord length from delta_c
% delta_c = 45;
delta_c = aero.delta_c0; %initial half-cone angle
chord = aero.rc0/tan(delta_c*pi/180);

max_alt = 150e3;    % maximum allowable altitude (m)
min_alt = 0;        % minimum allowable altitude (m)


% Base Input Structure
in0 = cvdma_in(max_alt,'venus',3); %create default cvdma at venus

% Get atmospheric profile
temp = load('./data/atm_data/atm_venus_mc.mat');
mc_atm = temp.mc;
in0.p.atm.table = temp.nom_table;
% copy nominal table to guidance system
in0.v.gnc.g.p_cvdma.planet.atm_table = in0.p.atm.table;

% Define Entry State
in0.s.traj.gamma_pp = x0.efpa*pi/180;
in0.s.traj.lat = x0.lat0*pi/180; % initial latitude, rad
in0.s.traj.lon = x0.lon0*pi/180; % initial longitude, rad
in0.s.traj.az = x0.az0*pi/180; % initial azimuth, deg->rad
in0.s.traj.vel_pp_mag = x0.v0*1e3;
in0.s.traj.alt = x0.h0 * 1e3;
in0.s.traj.t_ini = 0;
in0.s.traj.t_max = 1000;

% Convert scalars to inertial position and velocity vectors
[in0.s.traj.r_pci_ini, in0.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in0.s.traj.lat, in0.s.traj.lon, in0.s.traj.alt, ...
    in0.s.traj.gamma_pp, in0.s.traj.az, in0.s.traj.vel_pp_mag, ...
    in0.p.omega, in0.s.traj.t_ini, 0, in0.p.r_e, in0.p.r_p);

% Vehicle properties
in0.v.mp.m_ini = aero.m; % Initial mass (kg)
in0.v.mp.m_jettison = 0; % Jettisoned mass (kg) - unused for this guidance algorithm
in0.v.aero.area_ref = area_ini; % initial area difference before and after jettison (m^2)
in0.v.aero.mode = uint8(1); % constant lift and drag coefficients
in0.v.aero.cd = aero.cd0;
in0.v.aero.cl = 0;
in0.v.aero.nose_radius = aero.rn; % Nose radius (m)

% Trajectory Settings
in0.s.traj.rate = sim.traj_rate;           % double, Hz, integration rate
in0.v.gnc.g.p.rate = gnc.guid_rate;         % nd, guidance call rate (Hz)
in0.v.gnc.c.p.rate = in0.s.traj.rate;    % control rate = traj rate
in0.s.data_rate = sim.data_rate;
in0.v.gnc.n.p.rate = in0.s.traj.rate;   % navigation rate = traj rate

% Termination conditions
in0.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(1) = min_alt; % double, m, terminate at 0m altitude
in0.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

in0.s.term.or.type(2) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(2) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(2) = max_alt; % double, m, terminate at top of atmosphere
in0.s.term.or.direction(2) = int8(1); % int8, nd, positive crossing

% overall guidance settings
in0.v.gnc.g.p.dm_mode = uint8(5); % continuous drag modulation
in0.v.gnc.g.p.K_gain = 0.2;
in0.v.gnc.g.p.K_bounds = [0.1 2];
in0.v.gnc.g.p.Kflag = uint8(0); % no density correction on nominal

% CVDMA guidance settings
in0.v.gnc.g.p_cvdma.mode = uint8(1);    % npc
in0.v.gnc.g.p_cvdma.iter_max = uint8(1);
in0.v.gnc.g.p_cvdma.A_sens_atm = 0.25;
in0.v.gnc.g.p_cvdma.ha_tgt = gnc.ha_tgt*1e3;    %km -> m
in0.v.gnc.g.p_cvdma.ha_tol = 10e3;    %km -> m
in0.v.gnc.g.p_cvdma.cd_ini = aero.cd0;
in0.v.gnc.g.p_cvdma.rn = aero.rn;
in0.v.gnc.g.p_cvdma.delta_c_ini = delta_c*pi/180;
in0.v.gnc.g.p_cvdma.rc_ini = aero.rc0;
in0.v.gnc.g.p_cvdma.area_rate = gnc.area_rate;
in0.v.gnc.g.p_cvdma.rc_bounds = gnc.rc_bounds';
in0.v.gnc.g.p_cvdma.mass = aero.m;
in0.v.gnc.g.p_cvdma.c = chord;
in0.v.gnc.g.p_cvdma.t_max = in0.s.traj.t_max;
in0.v.gnc.g.p_cvdma.traj_rate = in0.s.traj.rate;

in0.v.gnc.g.p_cvdma.planet.alt_min = 55e3; % m
in0.v.gnc.g.p_cvdma.planet.alt_max = max_alt; % m

if in0.v.gnc.g.p_cvdma.rc_bounds(1) > in0.v.gnc.g.p_cvdma.rc_bounds(2)
    in0.v.gnc.g.p_cvdma.rc_bounds = flip(in0.v.gnc.g.p_cvdma.rc_bounds);
end

% Run nominal Trajectory
out = main1_mex(in0);
% out = main1(in0);         % debug mode
% 
% Compute parameters
if isnan(out.traj.alt(end))
    idxend = find(isnan(out.traj.alt),1)-1;
else
    idxend = length(out.traj.alt);
end

out.beta0 = aero.m/(aero.cd0*pi*aero.rc0^2);

rv = [out.traj.pos_ii(idxend,:), out.traj.vel_ii(idxend,:)]';
oe = rv2oe(rv,in0.p.mu);
ra = oe(1)*(1+oe(2));
haf = (ra - in0.p.r_e)/1000;

dv = out.traj.vel_ii_mag(1) - out.traj.vel_ii_mag(idxend);

out.haf = round(haf,5);                   %km
out.haf_err = round((haf - gnc.ha_tgt),5);    %km
out.dv = round(dv,5);
out.idxend = idxend;


% check monte carlo
if (mc.flag)
    out_nom = out;
    clear out;  % delete nominal trajectory
    out.nom = out_nom;  % nominal trajectory
    
%     in0.v.gnc.g.p.Kflag = uint8(mc.flag);   % turn on density corrector   
    % density corrector values way too high! turned off
    
    % preallocate vars
    len = in0.s.traj.t_ini:(1/in0.s.data_rate):in0.s.traj.t_max;
    len = length(len)+1;
    vmag = nan(len,mc.N);
    alt = vmag; fpa = vmag; t = vmag; 
    rho = vmag;
    ha = vmag; ha_err = vmag; 
    rc = vmag;
    K_dens = vmag;
    haf = nan(mc.N,1);
    haf_err = haf; 
    dv = haf;
%     tjett = nan;
    N = mc.N;
    sigs = mc.sigs;
    % pre-determine random atmospheres to reduce broadcast overhead
    temp = nan(1000, 7, N);
    for i = 1:N
        rnd = randi(1000);
        temp(:,:,i) = mc_atm(rnd).table;
    end
%     tic
    parfor i = 1:N
    % for i = 1:N
        in_mc = in0;        % copy nominal input struct
%         if mod(i,N/10) == 0
%             fprintf('Monte Carlo %i percent done\n',i*100/N);
%         end
        % apply disperse atmosphere
        
        in_mc.p.atm.table = temp(:,:,i);    % new atmospheric table
        
%         figure(); hold on
%         plot(in0.p.atm.table(:,2),in0.p.atm.table(:,1)./1000)
%         plot(in_mc.p.atm.table(:,2),in_mc.p.atm.table(:,1)./1000)
%         ylim([90 150])
        
        % apply input perturbations (only to traj. guid struct has nominal)
        in_mc = apply_dispersions(in_mc,sigs);
        
        % Convert scalars to inertial position and velocity vectors
        [in_mc.s.traj.r_pci_ini, in_mc.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
            in_mc.s.traj.lat, in_mc.s.traj.lon, in_mc.s.traj.alt, ...
            in_mc.s.traj.gamma_pp, in_mc.s.traj.az, in_mc.s.traj.vel_pp_mag, ...
            in_mc.p.omega, in_mc.s.traj.t_ini, 0, in_mc.p.r_e, in_mc.p.r_p);
        
%         out_mc = main1(in_mc);    % debugging purposes
        out_mc = main1_mex(in_mc);  % output struct for disperse case
        
        vmag(:,i) = out_mc.traj.vel_pp_mag./1000;   %km/s
        alt(:,i) = out_mc.traj.alt./1000;       %km
        fpa(:,i) = out_mc.traj.gamma_pp.*180/pi;    %deg
        t(:,i) = out_mc.traj.time;                  %s
        rho(:,i) = out_mc.traj.rho;     % kg/m3
        
        rc(:,i) = out_mc.g.cvdma.rc;
        
        ha(:,i) = out_mc.g.cvdma.ra./1000;
        ha_err(:,i) = out_mc.g.cvdma.delta_ra./1000;
        K_dens(:,i) = out_mc.g.K_dens;
        
        % get final index
        if isnan(out_mc.traj.alt(end))
            idxend = find(isnan(out_mc.traj.alt),1)-1;
        else
            idxend = length(out_mc.traj.alt);
        end
        
        % calculate apoapsis
        rv = [out_mc.traj.pos_ii(idxend,:), out_mc.traj.vel_ii(idxend,:)]';
        oe = rv2oe(rv,in0.p.mu);
        ra = oe(1)*(1+oe(2));
        err = (ra - in0.p.r_e)/1000;
        
        haf(i) = round(err,5);       %km, final apoapsis
        haf_err(i) = round((err - gnc.ha_tgt),5);    %km, apoapsis error
        dv(i) = round( out_mc.traj.vel_ii_mag(1) - out_mc.traj.vel_ii_mag(idxend), 5 );
    end % monte carlo loop
%     toc
    % store all monte carlo results into out struct
    out.traj.vmag = vmag;
    out.traj.alt = alt;
    out.traj.fpa = fpa;
    out.traj.t = t;
    out.traj.rho = rho;
    out.traj.rc = rc;
    out.g.ha = ha;
    out.g.ha_err = ha_err;
    out.g.K_dens = K_dens;
    out.haf = haf;
    out.haf_err = haf_err;
    out.dv = dv;
    
    
end % check mc flag


end % cv_dma_mars
