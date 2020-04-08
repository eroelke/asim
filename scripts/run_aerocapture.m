% run_aerocapture.m
% Aerocapture Master Function
%       Created Oct 2017, Evan Roelke    
%       Function to run any aerocapture script, any guidance mode
% 
%       *Optimized, Feb 2019, E. Roelke
%           ** unfinished
% 
% Inputs:
%     Simulation Structures
%       s:              simulation settings
%           s.p:            planet, string or number)
%           s.atm_mode      atmospheric mode, or specific table
%           s.traj_rate:    trajectory integration rate, Hz
%           s.data_rate:    data storage rate, Hz
%           s.t0:           initial time, s
%           s.tf:           final time, s
% 
%       x0:             entry conditions
%           x0.efpa:        entry flight path angle, deg
%           x0.v0:          entry (pcpf) velocity, km/s
%           x0.h_max:       max altitude, km
%           x0.h_min:       min altitude, km
%           x0.lat0:        initial latitude, deg
%           x0.lon0:        initial longitude, deg
%           x0.az0:         initial azimuth, deg
%       veh:            vehicle data
%           veh.rn:         nose radius, m
%           veh.rc:         backshell radii, m
%           veh.delta_c:    half-cone angle, deg
%           veh.m:          vehicle masses, kg
%           veh.cds:        drag coefficients
%           veh.cls:        lift coefficients
%           veh.aoa0:       initial angle of attack
%           veh.bank0:      initial bank angle
%           veh.aero_mode:  aerodynamics mode
%       gnc:            guidance, navigation, and control inputs
%           gnc.guid_rate:      guidance call rate, Hz
%           gnc.ctrl_rate:      control call rate, Hz
%           gnc.nav_rate:       navigation call rate, Hz
%           gnc.dm:             drag-modulation mode, uint8
%           gnc.prop_mode:      other guid mode (lift modulation), uint8
%           gnc.ha_tgt:         target apoapsis altitude, km
%           gnc.ha_tol:         apoapsis tolerance, km
%           
%       mc:             monte carlo struct
%           mc.N:                number of runs
%           mc.flag:            run monte carlo? (1)
%           mc.sigs:            uncertainty on entry conditions
% 
% 
% Outputs:
%     out_dat: output data structure
%         dv_pp: delta velocity, planet-relative
%         dv_ii: delta velocity, inertial-relative
%         fpa_corr: final flight path angle corridor
%         traj: relevant trajectory data for each simulation
%             alt: altitude data
%             time: temporal data
%             vel: velocity data
%             fpa: flight path angle
%         oe: orbital elements for each simulation
%             oe(1) = semi-major axis (sma)
%             oe(2) = eccentricity magnitude (ecc)
%             oe(3) = inclination (inc)
%             oe(4) = argument of periapsis (apr)
%             oe(5) = longitude of ascending node (lan)
%             oe(6) = true anomaly angle (taa)
%         haf: final apoapsis altitude (km)
%         
%         

function [out, dat_steep, dat_shallow, beta_ratio] = run_aerocapture(...
    sim_dat, planet, atm_mode, ...
    rcs, cDs, cLs, aoas, deltas, rns, m, aero_mode,  ...    %spacecraft geometry
    x_atm, sigs, ...                   %entry state vector and uncertanties
    fpa_corr_guess, ha_tgt, ha_tol, fpa_tol, ...
    guid_mode, guid_param, guid_trig, trig_tol)

addpath(genpath('C:\Users\Evan Roelke\Documents\research\asim'));

%define output structures
% out = define_out(sim_dat);   %see bottom of func

%% Planetary Data
bodies = {'venus','earth','mars','jupiter','saturn','titan','uranus','neptune'};
body_ind = find(strcmpi(bodies,planet) == 1);   %index of planetary bodies

atm_table = load('C:\Users\Evan\Documents\Academics\Research\asim\data\atm_venus_mc.mat');

% [in0, atm_table] = get_input_struct(body_ind,guid_mode);
in0 = sej_a_venus_in;
in0.p.atm.mode = uint8(atm_mode);        %planetary atmosphere mode

%% Entry State Vector
in0.s.traj.gamma_pp = x_atm(3)*pi/180;
in0.s.traj.alt = x_atm(1); % initial altitude, m
in0.s.traj.lat = x_atm(4)*pi/180; % initial latitude, rad
in0.s.traj.lon = x_atm(5)*pi/180; % initial longitude, rad
in0.s.traj.az = x_atm(6)*pi/180; % initial azimuth, deg->rad
in0.s.traj.vel_pp_mag = find_vel_pp_mag(x_atm(2), in0.s.traj.gamma_pp, in0.p.omega, (in0.p.r_e + in0.s.traj.alt)); % initial velocity, m/s

%% Convert scalars to inertial position and velocity vectors
[in0.s.traj.r_pci_ini, in0.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in0.s.traj.lat, in0.s.traj.lon, in0.s.traj.alt, ...
    in0.s.traj.gamma_pp, in0.s.traj.az, in0.s.traj.vel_pp_mag, ...
    in0.p.omega, in0.s.traj.t_ini, 0, in0.p.r_e, in0.p.r_p);

%% Simulation Settings
if sim_dat(3) ~= 0
    in0.s.traj.t_max = t_max; % double, s, maximum allowable time
end
in0.s.traj.rate = 100; % double, Hz, integration rate
in0.s.traj.t_ini = 0; % double, s, initial time
in0.s.data_rate = 1;

%% Termination conditions
in0.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(1) = int8(0); % double, m, terminate at 0m altitude
in0.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

in0.s.term.or.type(2) = int8(0); % int8, nd, crossing/equality condition
in0.s.term.or.var(2) = uint8(1); % uint8, nd, termination mode (altitude)
in0.s.term.or.value(2) = x_atm(1); % double, m, terminate at atmospheric exit
in0.s.term.or.direction(2) = int8(1); % int8, nd, positive crossing

%% Initial Vehicle data
in0.v.mp.m_ini = m(1);  % initial mass, kg
in0.v.aero.aoa = aoas(1)*pi/180;    %entry angle of attack, rad
in0.v.aero.cl = cLs(1); % entry lift coefficient
in0.v.aero.cd = cDs(1); % entry drag coefficient
in0.v.aero.area_ref = pi*rcs(1)^2; % entry ref area (m^2)
in0.v.aero.nose_radius = rns(1); % entry nose radius (m)
delta_c = deltas(1);      %entry spherecone half angle
beta(1) = in0.v.mp.m_ini/( in0.v.aero.cd * in0.v.aero.area_ref );

%% Aerocapture Settings
fpa_tol = fpa_tol*pi/180;   %degrees to radians
fpa_corr_guess = fpa_corr_guess.*pi/180;


%% get aerodynamics data
switch aero_mode
    case 0  %constant aerodynamics
        in0.v.aero.mode = uint8(1); % aerodynamics mode = constant cl and cd
    case 1  %mach-dependent table lookup
        in0.v.aero.mode = uint8(2); % aerodynamics mode = table look up
        
        if delta_c == 45
            temp = load(['.' filesep 'data' filesep 'aero_sphere_cone_45deg.mat']);
        elseif delta_c == 60
            temp = load(['.' filesep 'data' filesep 'aero_sphere_cone_60deg.mat']);
        elseif delta_c == 70
            temp = load(['.' filesep 'data' filesep 'aero_sphere_cone_70deg.mat']);
        else
            fprintf('Invalid aerodynamics mode input...\nLoading 70deg SphereCone default\n');
            temp = load(['.' filesep 'data' filesep 'aero_sphere_cone_70deg.mat']);
        end
        in0.v.aero.table = temp.table;
    case 2  %mach and aoa dependent table lookup
        
        %%TO DO
        
    otherwise
        error('Incorrect Aerodynamics Mode Input');
end % aero_mode

%% run nominal trajectory
% in_nom = in0;
% in_nom.p.atm.table = atm_table.nom_table;
% 
% switch guid_mode
%     
%     case 0  %no guidance
%         % TO DO
%     case 1  %manual SEJ
%         % TO DO        
%     case 2  %auto SEJ
%         in_nom.v.mp.m_jettison = m(1) - m(2);
%         in_nom = set_sej_guidance(in_nom,x_atm,rcs);
%         
%         if fpa_corr_guess ~= 0
%             in_temp = in_nom;
%             in_temp.v.gnc.g.p.prop_mode = uint8(1); %no guidance
%             in_temp.v.aero.area_ref = pi*rcs(1)^2;
%             [fpa_corr(1), ~, dat_shallow, ~, ~] = ... 
%                 find_fpa_given_ha(fpa_corr_guess, ha_tgt, in_temp, fpa_tol, ha_tol);
%             
%             in_temp.v.aero.area_ref = pi*rcs(2)^2;
%             in_temp.v.mp.m_ini = m(2);
%             [fpa_corr(2), ~, dat_steep, ~, ~] = ... 
%                 find_fpa_given_ha(fpa_corr_guess, ha_tgt, in_temp, fpa_tol, ha_tol);
%             fpa_corr = fpa_corr.*180/pi;
%             if abs(x_atm(3)) < abs(fpa_corr(1)) || abs(x_atm(3)) > abs(fpa_corr(2))
%                 fprintf('WARNING: Entry FPA outside of corridor!\n')
%                 
%                 ans = str(input('Do you want to adjust FPA_atm? (y/n): '));
%                 if ( strcmp(ans,'y') )
%                    x_atm(3) = mean(fpa_corr);
%                 else
%                     return;
%                 end
%                 
%             end
%         else
%             fpa_corr = nan(1,2);
%             dat_shallow = nan;
%             dat_steep = nan;
%         end
%         
%         beta(2) = m(2)/( cDs(2)*pi*rcs(2)^2 );
%         beta_ratio = beta(2)/beta(1);
%         
%         % run trajectory
%         in_nom.v.gnc.g.p_sej_a.tgt_ap = ha_tgt; %target apoapsis, m
%         dat_nom = main1_mex(in_nom);
%         
%         %check capture
%         out.nom.orbit = get_ac_oe(dat_nom,in_nom.p.mu,in_nom.p.r_e);
%         haf_corr = [ha_tgt - ha_tol, ha_tgt + ha_tol]; %tolerable ha
%         if out.nom.orbit.haf < haf_corr(1) || out.nom.orbit.haf > haf_corr(2)
%             fprintf('Unsuccessful Capture. haf = %4.2f m, diff = %4.2f m\n', ... 
%                 out.nom.orbit.haf, abs(out.nom.orbit.haf - ha_tgt) );
% %             out.nom.orbit.haf = nan;
% %             return;
%         else
%             fprintf('Successful Capture! haf = %4.2f m\n',out.nom.orbit.haf);
%         end
%         
%         %get output data
%         out.nom.traj = dat_nom.traj;
%         
%         
%         
%     case 4  %manual continous drag
%         % TO DO
%         N = length(rcs);    %number of control segments
%         beta1 = beta(1);    beta = nan(N,1);
%         beta_ratio_rel = beta;
%         beta_ratio_abs = beta;
%         beta(1) = beta1;
%         for i = 1:N
%            
%             beta(i+1) = m(i+1)/( cd(i+1)*pi*rc(i+1)^2 );
%             beta_ratio_rel(i) = beta(i+1)/beta(i);  %relative to last beta
%             beta_ratio_abs(i) = beta(i+1)/beta(1);  %relative to initial beta
%             
%         end
%         beta_ratio = [beta_ratio_rel', beta_ratio_abs'];
%         
%     case 5  %manual continuous lift modulation
%         % TO DO
%     case 6  %auto continious lift modulation
%         % TO DO
%     otherwise
%         error('Incorrect Guidance Mode')
% end


% out_nom = [];

%% run monte carlo trajectory
N = sim_dat(2);  % number of monte carlo cases
if N > 1    %
% in_mc.p.atm.table = atm_table.mc;   %load 1000 different tables
% fprintf('Starting Monte Carlo Simulation...\n')

% preallocate
% out.traj.alt = nan(in0.s.traj.t_max,N);
% out.traj.vel = nan(in0.s.traj.t_max,N);
% out.traj.time = nan(in0.s.traj.t_max,N);
% out.traj.heat = nan(in0.s.traj.t_max,N);
% out.traj.q = nan(in0.s.traj.t_max,N);
% out.traj.g = nan(in0.s.traj.t_max,N);
% out.traj.fpa = nan(in0.s.traj.t_max,N);

out.orbit.haf = nan(N,1);
out.orbit.dv = nan(N,1);    % aerocapture velocity benefit

for i = 1:N
    if mod(i,5) == 0
        percent = (i/N) * 100;
        fprintf('Run %i/%i. %4.2f percent done...\n',i,N,percent);
    end
    in_mc = in_nom;
    temp = randi(1000);
    in_mc.p.atm.table = atm_table.mc(temp).table;  % load random atmosphere
    
    %apply dispersions to entry state, etc
    in_mc = apply_dispersions(in_mc,sigs);
    
    % run simulation
    dat_mc = main1_mex(in_mc);
    
    out.mc(i).traj = dat_mc.traj;
    
    out.mc.traj(i).alt = dat_mc.traj.alt;
    out.mc.traj(i).vel = dat_mc.traj.vel_pp_mag;
    out.mc.traj(i).time = dat_mc.traj.time;
    out.mc.traj(i).q = dat_mc.traj.q;
    out.mc.traj(i).g = dat_mc.traj.g_loading;
    out.mc.traj(i).fpa = dat_mc.traj.gamma_pp;
    out.mc.traj(i).heat_rate = dat_mc.traj.heat_rate;
        
%     out.traj.alt(:,i) = dat_mc.traj.alt;  % altitude
%     out.traj.vel(:,i) = dat_mc.traj.vel_pp_mag; % planet relative velocity
%     out.traj.time(:,i) = dat_mc.traj.time;  % simulation time, sec
%     out.traj.heat(:,i) = dat_mc.traj.heat_rate;   % peak heat rate
%     out.traj.q(:,i) = dat_mc.traj.q;
%     out.traj.g(:,i) = dat_mc.traj.g_loading;
%     out.traj.fpa(:,i) = dat_mc.traj.gamma_pp;
%     
    [out.orbit.haf(i),out.orbit.dv(i)] = ... 
        get_ac_result(dat_mc,in_mc.p.mu,in_mc.p.r_e);
    
    
    %check capture
%     out.mc(i).orbit = get_ac_oe(dat_mc,in_mc.p.mu,in_mc.p.r_e);
    
%     if out.mc(i).orbit.haf < haf_corr(1) || out.mc(i).orbit.haf > haf_corr(2)
%         fprintf('Unsuccessful Capture. haf = %4.2f m, diff = %4.2f m\n', ...
%             out.mc(i).orbit.haf, abs(out.mc(i).orbit.haf - ha_tgt) );
%         out.mc.orbit(i,:).haf = nan;
%         return;
%     else
%         fprintf('Successful Capture! haf = %4.2f m\n', ...
%             out.mc(i).orbit.haf );
    
%     end
    
    
end % N_mc loop


end % check monte carlo

%% final output settings
out.fpa_corr = fpa_corr;


end %run_aerocapture


function in = set_sej_guidance(in,x_atm,rcs)
% set single event jettison guidance parameters

in.v.gnc.g.p_sej_a.alt_min = 0; % m
in.v.gnc.g.p_sej_a.alt_max = x_atm(1); % m
in.v.gnc.g.p_sej_a.cd = in.v.aero.cd; % nd
in.v.gnc.g.p_sej_a.cl = in.v.aero.cl; % nd
in.v.gnc.g.p_sej_a.mass = in.v.mp.m_ini; % kg
in.v.gnc.g.p_sej_a.mass_jettison = in.v.mp.m_jettison;
in.v.gnc.g.p_sej_a.area_ref = pi.*[rcs(1)^2 rcs(2)^2 0]';

in.v.gnc.g.p_sej_a.atm_table = in.p.atm.table;

end


function [out] = define_out(sim_dat)

orbit = struct(...
    'rv', nan, ...
    'oe',nan, ...
    'ra',nan, ...
    'haf',nan);

nom = struct(...
    'traj', nan, ...
    'fpa_corr', nan, ...
    'orbit', orbit);

for i = 1:sim_dat(2)
    mc(i) = nom;
end

out = struct(...
    'nom', nom, ...
    'mc', mc, ...
    'fpa_corr', nan(1,2));

end

function [oe] = get_ac_oe(dat,mu,r_e)

%get final trajectory index
if ( isnan(dat.traj.alt(end)) )
    idxend = find(isnan(dat.traj.alt),1)-1;
else
    idxend = length(dat.traj.alt);
end

rv = [dat.traj.pos_ii(idxend,:), dat.traj.vel_ii(idxend,:)]';
oe = rv2oe(rv,mu);
ra = oe(1)*(1+oe(2));   %m
haf = (ra - r_e);       %m

dv = dat.traj.vel_ii_mag(1) - dat.traj.vel_ii_mag(idxend);


oe = struct(...
    'rv', rv, ...
    'oe', oe, ...
    'ra', ra, ...
    'haf', haf, ...
    'dv', dv);

end



function [haf,dv] = get_ac_result(dat,mu,r_e)

%get final trajectory index
if ( isnan(dat.traj.alt(end)) )
    idxend = find(isnan(dat.traj.alt),1)-1;
else
    idxend = length(dat.traj.alt);
end

rv = [dat.traj.pos_ii(idxend,:), dat.traj.vel_ii(idxend,:)]';
oe = rv2oe(rv,mu);
ra = oe(1)*(1+oe(2));   %m
haf = (ra - r_e);       %m

dv = dat.traj.vel_ii_mag(1) - dat.traj.vel_ii_mag(idxend);


oe = struct(...
    'rv', rv, ...
    'oe', oe, ...
    'ra', ra, ...
    'haf', haf, ...
    'dv', dv);

end


