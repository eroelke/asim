%Titan Traj
%
%   Inputs:
%       winds: sinusoidal ('sin') or regular ('reg') wind profiles
%       run_disp: run disperse trajectories, (0 or 1)
%       vehicle: vehicle configuration - 'sc','5B','HB','SO2P','IXV','AOTV'
%           'sc' - spherecone (can include HIAD design)
%           '5B' - AMOOS 5B
%           'HB' - AMOOS HB
%           'SO2P' - Saturn Dual Orbiter Vehicle
%           'IXV' - ESA IXV
%           'AOTV' - Inflatable Chine AOTV
%       n: number of disperse trajectories to run
%       t_max: maximum time limit on simulation
%       v_atm: atmospheric entry velocity, m/s
%       m_atm: entry mass - scales the reference area for aerodynamic forces, kg
%       fpa_atm: flight path angle at atmospheric interface, deg
%       
%       rn: nose radius, for spherecone vehicle ('sc'), m
%       rc: reference area radius for spherecone vehicle (can include HIAD config), m
%       aoa: angle of attack, determines L/D ratio, deg
%       theta: sphere cone half angle, deg
%       cl: hypersonic lift coefficient
%       cd: hypersonic drag coefficient
%
%       h_sig: atmospheric interface 1sigma uncertainty, m (Huygens 30km)
%       v_sig: velocity at atmospheric interface 1sigma uncertainty, m/s (Huygens 10.15 m/s 3sigma)
%       fpa_sig: flight path angle at atmospheric interface 1sigma uncertainty, deg (Huygens 0.85 degrees)
%
%       num_chutes: 1 or 2 parachute system
%       d_chutes: diameter of parachute(s), (m)
%       v_chute_trigger: velocity at which to deploy the parachute, (m/s)
%
%   Outputs:
%       nom_v_sd: nominal impact velocity (m/s)
%       nom_dr: nominal downrange (km)
%       nom_cr: nominal crossrange (km)
%       nom_peak_heat: nominal peak heat rate (W/cm2)
%       nom_heat_load: nominal heat load (J/cm2)
%       nom_peak_g: nominal peak deceleration in Earth Gs
%       nom_time: nominal descent time (hrs)
%       max_sd: maximum impact velocity for disperse traj, m/s
%       min_sd: minimum impact velocity for disperse traj, m/s
%       med_sd: median impact velocity, m/s
%       avg_sd: average impact velocity, m/s
%       std_sd: standard deviation impact velocity, m/s
%       mc_time: maximum descent time in disperse traj, hrs
%       dr_range: downrange footprint, km
%       cr_range: crossrange, footprint, km
%       mc_peak_heat: maximum peak heat rate, W/cm2
%       mc_heat_load: maximum heat load, J/cm2
%       mc_peak_g: maximum peak deceleration, Earth G's
%
% Created By: Evan Roelke OCT 2016
% 

function [nom_v_sd, nom_dr, nom_cr, nom_peak_heat, nom_heat_load, nom_peak_g, nom_time, ...
    max_sd,min_sd,med_sd,avg_sd,std_sd, ...
    mc_time,dr_range,cr_range,mc_peak_heat,mc_heat_load,mc_peak_g] = TitanTraj ... 
    (winds,run_disp,vehicle, n, t_max, ...    %wind profile, disperse, vehicle configuration
    v_atm,m_atm,fpa_atm, ...        %nominal entry state vector
    rn,rc,aoa,theta,cl,cd, ...            %spherecone geometry (for vehicle 'sc'
    h_sig,v_sig,fpa_sig, ...     %uncertanties for error propogation    
    num_chutes, d_chutes, v_chute_trigger)           %parachute(s)
%#codegen
format longg;
%% Define input structure
    in = default_in;
%% Planetary data
    in.p = get_planetary_data(6,3); % Load settings for Titan (6), model with winds (3)
    if strcmp(winds,'sin') == 1
        atm_tables = load([ '.' filesep 'data' filesep 'atm_titan_gram_10000mc_sin_winds.mat']);
        in.p.atm.table = atm_tables.atm_mc_sin_tables.nom_table;
    else
        atm_tables = load([ '.' filesep 'data' filesep 'atm_titan_gram_10000mc_reg_winds.mat']);
        in.p.atm.table = atm_tables.atm_mc_reg_tables.nom_table; winds = 'regular';
    end


%% Simulation data
% Termination conditions
    in.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
    in.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
    in.s.term.or.value(1) = 0; % double, m, terminate at 0 km altitude (sea level)
    in.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

% Integration parameters
    in.s.traj.rate = 10; % double, Hz, integration rate
    in.s.traj.t_ini = 0; % double, s, initial time
    in.s.traj.t_max = t_max; % double, s, maximum allowable time
    in.s.data_rate = 1;

% Initial state (Eastbound equatorial)
    in.s.traj.alt = 1270e3; % atmospheric interface altitude, m
    in.s.traj.vel_pp_mag = v_atm; % initial velocity, m/s
    
    in.s.traj.lat = 81 * pi/180; % initial latitude (72.4deg north), rad
    in.s.traj.lon = -338*pi/180; % initial longitude (338deg west), rad
    in.s.traj.gamma_pp = fpa_atm*pi/180; % initial flight-path, rad
    in.s.traj.az = 88.89*pi/180; % initial azimuth, rad, Huygens Probe
%     in.s.traj.az = 0;   %initial azimuth, rad

% Convert scalars to inertial position and velocity vectors
    [in.s.traj.r_pci_ini, in.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
        in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
        in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
        in.p.omega, in.s.traj.t_ini, 0, in.p.r_e, in.p.r_p);

%% Vehicle data
if strcmp(vehicle,'sc') == 1
    in.v.mp.m_ini = m_atm; % kg, vehicle mass at entry
%     if aoa == 0
%         if theta == 70
%             in.v.aero.mode = uint8(2); % nd, table look up
%             temp = load( ['.' filesep 'data' filesep 'aero_sphere_cone_70deg.mat'] );
%             in.v.aero.table = temp.table;
%         elseif theta == 60
%             in.v.aero.mode = uint8(2); % nd, table look up
%             temp = load( ['.' filesep 'data' filesep 'aero_sphere_cone_60deg.mat'] );
%             in.v.aero.table = temp.table;
%         elseif theta == 45
%             in.v.aero.mode = uint8(2); % nd, table look up
%             temp = load( ['.' filesep 'data' filesep 'aero_sphere_cone_45deg.mat'] );
%             in.v.aero.table = temp.table;
%         else
%             in.v.aero.mode = uint8(1);  %user defined cL and cD
%         end
%     else
        in.v.aero.mode = uint8(1);  %constant cl and cd
        in.v.aero.cl = cl;
        in.v.aero.cd = cd;
%     end
    
    in.v.aero.nose_radius = rn;  %nose radius, meters
    in.v.aero.area_ref = pi*(rc)^2;  %reference area
    in.v.aero.aoa = aoa; %angle of attack, radians
    
    in.v.decel.K = 0; % no dispersions (nominal)
    in.v.decel.mode = uint8(1); %user-input
    in.v.decel.d = d_chutes(1); % first parachute diameter, m
    in.v.decel.cd = 0.55;
    
elseif strcmp(vehicle,'5B') == 1
    %AMOOS 5B
    in.v.mp.m_ini = 2800;
    in.v.aero.mode = uint8(1); % nd, aerodynamics mode
        in.v.aero.area_ref = pi*(2.2352)^2; % m^2, aerodynamic reference area based 
            %on AMOOS reresentative packaging (18.3m long) original idea of X-37B is 3m high 
        in.v.aero.nose_radius = 0.88; % m, nose radius: AMOOS packaging (for 18.3 long)
        in.v.aero.cd = 3.285; % nd, cD at cL_max
        in.v.aero.cl = 1.979; % nd, cL_max
        in.v.aero.aoa = 45 * pi/180;    %trim angle of attack at cL_max
        in.s.traj.aoa_ini = in.v.aero.aoa;
        
        in.v.decel.d = 7;
elseif strcmp(vehicle,'HB') == 1
    %AMOOS HB
    in.v.mp.m_ini = 2800;
    scale = 0.45;       %size scale factor
    in.v.aero.mode = uint8(1); % nd, aerodynamics mode
    in.v.aero.area_ref = 14.864*scale^2; % m^2, aerodynamic reference area based
    %on AMOOS reresentative packaging (18.3m long) original idea of X-37B is 3m high
    in.v.aero.nose_radius = 3.7253*scale; % m, nose radius: AMOOS packaging (for 18.3 long)
    in.v.aero.cd = 3.874; % nd, cD at cL_max
    in.v.aero.cl = 2.003; % nd, cL_max
    in.v.aero.aoa = 45 * pi/180;    %trim angle of attack at cL_max
    
    in.v.decel.d = 7;
elseif strcmp(vehicle,'SO2P') == 1
    %SO2P
    in.v.decel.d = 7;
elseif strcmp(vehicle,'AOTV') == 1
    %inflatable Chine AOTV
    in.v.decel.d = 7;
elseif strcmp(vehicle,'IXV') == 1
    %ESA IXV
    in.v.mp.m_ini = 2800;
    in.v.aero.mode = uint8(1); % nd, aerodynamics mode
    in.v.aero.area_ref = 7.26; % m^2, aerodynamic reference area
    % in.v.aero.nose_radius = ; % m, nose radius:
    in.v.aero.cd = 1.515; % nd, cD at cL_max
    in.v.aero.cl = 1.0605; % nd, cL_max from published L/D value
    %         in.v.aero.cl = 0.06;    %cL_max from data
    in.v.aero.aoa = -25*pi/180;    %trim angle of attack at cL_max
    
    in.v.decel.d = 7;
else
    error('Vehicle Specification Incorrect')
end

%% GNC
% Navigation 
    in.v.gnc.n.p.mode = uint8(2); % Use Markov process/ECRV error model
    in.v.gnc.n.p.rate = 20; % Hz
    in.v.gnc.n.p.seed = uint32(1); % nd, Error model seed
    in.v.gnc.n.p.tau = 3600; % s, time constant
    in.v.gnc.n.p.omega = in.p.omega; % rad/s, planet angular velocity vector
    in.v.gnc.n.p.r_e = in.p.r_e; % m, planet equatorial radius
    in.v.gnc.n.p.r_p = in.p.r_p; % m, planet polar radius

% Guidance
    %Bank
    in.v.gnc.g.p.rate = 0.5; % Hz
    in.v.gnc.g.p.bank_mode = uint8(2); % nd, constant bank rate
    in.v.gnc.g.p.const_bank_rate = 4*pi/60; % rad/s, 2 RPM
    in.v.gnc.g.p.const_bank = in.v.gnc.g.p.const_bank_rate; 
    in.v.gnc.g.p.tgt_lat = 81*pi/180; % rad (81deg N - Kraken Mare center)
    in.v.gnc.g.p.tgt_lon = 110*pi/180;        %250 W = 110 E

% Control
    in.v.gnc.c.p.bank_mode = uint8(3);  %accel and rate-limit control
    in.v.gnc.c.p.rate = 20; % Hz
    in.v.gnc.c.p.decel_mode = uint8(1);
    in.v.gnc.c.s.aoa = in.v.aero.aoa;
    in.v.gnc.c.s.bank_rate = 4*pi/60;

    in.v.gnc.g.p.aoa_mode = uint8(1);
    in.v.gnc.g.p.const_aoa = in.v.aero.aoa; 
    in.s.traj.aoa_ini = in.v.aero.aoa;


%% Nominal Trajectory Analysis
    terminate = 0;
    if num_chutes == 2
        %deployment parameters
        mach_deploy = 1.5;
        v_deploy = v_chute_trigger;
        d_pf = d_chutes(2);
        d_pc = d_chutes(1);
        [~,~,~,terminate,~,~,~,vel,~,time,dr,cr,heat_rate,~,decel,~,~,~,~,~,~] = ...
            TitanTraj_TwoChutes(in,mach_deploy,v_deploy,d_pf,d_pc);
    elseif num_chutes == 1
        %Slender code
        mach_deploy = 1.5;
        d_chute = d_chutes(1);
        [~,~,~,~,vel,~, ...
            time,dr,cr,heat_rate,~,decel,~,~,~,~,~,~] = ... 
            TitanTraj_OneChute(in,mach_deploy,d_chute);
    else
        error('Only Single or Double Parachute Systems are considered')
    end
    
    if terminate == 1
        disp('Vehicle Reached Ground Before Second Chute Deployed')
    end
%nominal output values
    nom_v_sd = vel(end);
    nom_dr = dr(end)/1000;
    nom_cr = cr(end)/1000;
%     nom_lat = lat(end)*180/pi;
%     nom_lon = 360 - lon(end)*180/pi;
    nom_peak_heat = max(heat_rate)/(100^2);        %peak heating rate, W/cm2
    nom_peak_g = max(decel);           %peak deceleration, earth g's
%     nom_peak_q = max(q);                   %peak dynamic pressure, N/m2
    nom_time = time(end)/3600;     %descent time, hrs
    nom_heat_load = trapz(time,heat_rate)/(100^2);  %J/cm2
    
%% Disperse Trajectories
if run_disp == 1
    in.s.traj.run_disp = 1;
else
    keyboard
end
% n = 2000;   %number of trajectories
if in.s.traj.run_disp == 1;
    %initialize all variables
    %initialize output variables
    disp_vel = nan(in.s.traj.t_max + 2,n);    %velocity magnitude
    disp_alt = nan(in.s.traj.t_max + 2,n);    %altitude
    disp_time = nan(in.s.traj.t_max + 2,n);     %time
    disp_dr = nan(in.s.traj.t_max + 2,n);     %downrange
    disp_cr = nan(in.s.traj.t_max + 2,n);     %crossrange
    disp_heat_rate = nan(in.s.traj.t_max + 2,n);     %heat rate
    disp_q = nan(in.s.traj.t_max + 2,n);     %dynamic pressure
    disp_decel = nan(in.s.traj.t_max + 2,n);     %decleration g_loading
    disp_lat = nan(in.s.traj.t_max + 2,n);    %latitude
    disp_lon = nan(in.s.traj.t_max + 2,n);    %longitude
    disp_mass = nan(in.s.traj.t_max + 2,n);     %mass
    disp_drag = nan(in.s.traj.t_max + 2,n);  %total drag force (vehicle + chute)
    disp_az = nan(in.s.traj.t_max + 2,n);      %azimuth
    disp_gamma = nan(in.s.traj.t_max + 2,n);     %fpa
    disp_mach_ind = nan(1,n);                   %index when pilot chute deploys
    disp_ind2 = nan(1,n);                   %index when parafoil deploys
    disp_ind3 = nan(1,n);                       %index when vehicle reaches the ground
    
    %initailize input variables for random sampling
    %two parachute system
        disp_mach_pc = nan(1,n);
        disp_v_pf = nan(1,n);
        disp_terminate = nan(1,n);
    %single parachute system
        disp_mach_deploy = nan(1,n);
        
        %initialize output variables
        v_sd = nan(n,1);              %splashdown velocity
        peak_heat_mc = nan(n,1);        %peak heating rate, W/cm2
        heat_load_mc = nan(n,1);        %heat load in each case
        peak_g_mc = nan(n,1);           %peak deceleration, earth g's
        peak_q_mc = nan(n,1);                   %peak dynamic pressure, N/m2
        time_sd = nan(n,1);             %descent time each case, hrs
        dr_sd = nan(n,1);                  %downrange at splashdown. km
        cr_sd = nan(n,1);                  %crossrange at splashdown, km
    
    %other input variables
    tic;
    for i = 1:n
        i
        %apply randomly sampled uncertainties
        mc_val = randi(1000);     %generate random integer between 1 and 1000
        in_disp = apply_dispersions(in,h_sig,v_sig,fpa_sig);
        
        if strcmp(winds,'sin') == 1
            in_disp.p.atm.table = atm_tables.atm_mc_sin_tables.mc(mc_val).table;
        else
            in_disp.p.atm.table = atm_tables.atm_mc_reg_tables.mc(mc_val).table;
        end
        
        [in_disp.s.traj.r_pci_ini, in_disp.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
            in_disp.s.traj.lat, in_disp.s.traj.lon, in_disp.s.traj.alt, ...
            in_disp.s.traj.gamma_pp, in_disp.s.traj.az, in_disp.s.traj.vel_pp_mag, ...
            in_disp.p.omega, in_disp.s.traj.t_ini, 0, in_disp.p.r_e, in_disp.p.r_p);
        
        if num_chutes == 2
            disp_mach_pc(i) = mach_deploy + normrnd(0,0.01);    %mach number at pilot chute deployment +- 0.01 1 sigma
            disp_v_pf(i) = v_deploy + normrnd(0,1);             %velocity at parafoil deployment +- 1 m/s 1 sigma
            
            [~,~,~,disp_terminate(i),disp_mach_ind(i),disp_ind2(i),disp_ind3(i),disp_vel(:,i),disp_alt(:,i), ...
                disp_time(:,i),disp_dr(:,i),disp_cr(:,i),disp_heat_rate(:,i), disp_q(:,i),disp_decel(:,i), ...
                disp_lat(:,i),disp_lon(:,i),disp_mass(:,i),disp_drag(:,i),disp_az(:,i),disp_gamma(:,i)] = ...
                TitanTraj_TwoChutes(in_disp,disp_mach_pc(i),disp_v_pf(i),d_pf,d_pc,in.s.traj.run_disp);
        elseif num_chutes == 1
            disp_mach_deploy(i) = mach_deploy + normrnd(0,0.01);    %mach number at chute deployment +- 0.01 1 sigma
            
            [~,~,disp_mach_ind(i),disp_ind2(i),disp_vel(:,i),disp_alt(:,i), ...
                disp_time(:,i),disp_dr(:,i),disp_cr(:,i),disp_heat_rate(:,i), disp_q(:,i),disp_decel(:,i), ...
                disp_lat(:,i),disp_lon(:,i),disp_mass(:,i),disp_drag(:,i),disp_az(:,i),disp_gamma(:,i)] = ...
                TitanTraj_OneChute(in_disp,disp_mach_deploy(i),in.v.decel.d);
        else
            error('Only Single or Double Parachute Systems are considered')
        end
        
        %find end of trajectory
        sd_ind = find(isnan(disp_time(:,i)) == 1);
        sd_ind = sd_ind(1) - 1;
        
        v_sd(i) = disp_vel(sd_ind,i)';              %splashdown velocity
        peak_heat_mc(i) = max(disp_heat_rate(:,i))/(100^2)';        %peak heating rate, W/cm2
        heat_load_mc(i) = trapz(disp_time(1:sd_ind,i),disp_heat_rate(1:sd_ind,i)); %heat load in each case
        peak_g_mc(i) = max(disp_decel(:,i))';           %peak deceleration, earth g's
        peak_q_mc(i) = max(disp_q(:,i))';                   %peak dynamic pressure, N/m2
        time_sd(i) = disp_time(sd_ind,i)/3600';             %descent time each case, hrs
        dr_sd(i) = disp_dr(sd_ind,i)/1000';                  %downrange at splashdown. km
        cr_sd(i) = disp_cr(sd_ind,i)/1000';                  %crossrange at splashdown, km
%         lat_sd(i) = disp_lat(sd_ind,i)';                %latitude at splashdown, rad
%         lon_sd(i) = disp_lon(sd_ind,i)';                %longitude at splashdown, rad
        if disp_alt(sd_ind,i) > 1
            disp('vehicle did not reach the ground')
            keyboard
        end
        clear in_disp %clear data
    end
    toc;
    
    %calculate velocity stats
    max_sd = max(v_sd);       %maximum splashdown speed
    min_sd = min(v_sd);       %mimimum splashdown speed
    avg_sd = mean(v_sd);      %average splashdown speed
    med_sd = median(v_sd);      %median splashdown speed
    std_sd = std(v_sd);       %standard deviation splashdown speed
    
    %downrange stats
    max_dr = max(dr_sd);      %maximum downrange(km) from entry position
    min_dr = min(dr_sd);       %minimum downrange(km) from entry position
%     avg_dr = mean(dr_sd);      %mean downrange(km) from entry position
%     med_dr = median(dr_sd);     %median downrange(km) from entry
    dr_range = max_dr - min_dr;     %downrange landing ellipse
    
    %crossrange stats
    max_cr = max(cr_sd);      %maximum crossrange(km) from entry position
    min_cr = min(cr_sd);       %minimum crossrange(km) from entry position
%     avg_cr = mean(cr_sd);      %mean crossrange(km) from entry position
%     med_cr = median(cr_sd);     %median crossrange(km) from entry
    cr_range = max_cr - min_cr;     %crossrange landing ellipse
    
%     %latitude, longitude
%     lat_sd = lat_sd.*180/pi;            %degrees N
%     avg_lat = mean(lat_sd);
%     std_lat = std(lat_sd);
%     
%     lon_sd = 360 - lon_sd.*180/pi;      %degrees W
%     avg_lon = mean(lon_sd);
%     std_lon = std(lon_sd);
    
    %other
    mc_time = max(time_sd);
    mc_peak_heat = max(peak_heat_mc);
    mc_heat_load = max(heat_load_mc)/(100^2);
    mc_peak_g = max(peak_g_mc);
    
end