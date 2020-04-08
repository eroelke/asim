% Titan Submarine_6m
%   Create input structure for guided entry at Titan. Vehicle based
%   on idea of landing 6m long Titan Submarine in Kraken Mare,
%   subject to current launch vehicle maximum diameter constraints (Atlas
%   V)
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Major Revision History:
%   *Created MAR 2016 - E. ROELKE

%#codegen
clear; 
clc;
% clf;
format longg;
%% Define input structure
    in = default_in;

%% Planetary data
    in.p = get_planetary_data(6,3); % Load settings for Titan (6), model with winds (3)
    atm_tables = load([ '.' filesep 'data' filesep 'atm_titan_gram_10000mc_fixed_winds.mat']);
    in.p.atm.table = atm_tables.fixed_atm_mc_tables.nom_table;
    
%% Simulation data
% Termination conditions
    in.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
    in.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
    in.s.term.or.value(1) = 0; % double, m, terminate at 0 km altitude (sea level)
    in.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

% Integration parameters
    in.s.traj.rate = 10; % double, Hz, integration rate
    in.s.traj.t_ini = 0; % double, s, initial time
    in.s.traj.t_max = 20000; % double, s, maximum allowable time
    in.s.data_rate = 1;

% Initial state (Eastbound equatorial),
    in.s.traj.alt = 1270e3; % atmospheric interface altitude, m
    in.s.traj.lat = 73.4 * pi/180; % initial latitude (72.4deg north), rad
    % in.s.traj.lon = 22*pi/180; % initial longitude (3dr_ran38deg west), rad
    in.s.traj.lon = -338*pi/180; % initial longitude (338deg west), rad
    in.s.traj.vel_pp_mag = 6000; % initial velocity, m/s
    in.s.traj.gamma_pp = -65.4*pi/180; % initial flight-path, rad (Labreton ref)
    in.s.traj.az = 88.89*pi/180; % initial azimuth, rad
    tgt_lat = 80*pi/180; % rad (70deg N - Ligeia Mare center)
    % tgt_lon = -315*pi/180; % rad (315deg W - Kraken Mare center)?
    tgt_lon = -70*pi/180;        %70 W = 290E

% Convert scalars to inertial position and velocity vectors
    [in.s.traj.r_pci_ini, in.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
        in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
        in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
        in.p.omega, in.s.traj.t_ini, 0, in.p.r_e, in.p.r_p);

%% Vehicle data
% Mass properties
    in.v.mp.m_ini = 2800;   %kg, 6m long sub entry mass

% Aerodynamics
%     % AMOOS config 5B constant cL and cD
%         in.v.aero.mode = uint8(1); % nd, aerodynamics mode
%         in.v.aero.area_ref = pi*(2.2352)^2; % m^2, aerodynamic reference area based 
%             %on AMOOS reresentative packaging (18.3m long) original idea of X-37B is 3m high 
%         in.v.aero.nose_radius = 0.88; % m, nose radius: AMOOS packaging (for 18.3 long)
%         in.v.aero.cd = 3.285; % nd, cD at cL_max
%         in.v.aero.cl = 1.979; % nd, cL_max
%         in.v.aero.aoa = 45 * pi/180;    %trim angle of attack at cL_max
%         in.s.traj.aoa_ini = in.v.aero.aoa;

    %AMOOS config HB constant cL and cD
        scale = 0.45;       %size scale factor
        in.v.aero.mode = uint8(1); % nd, aerodynamics mode
        in.v.aero.area_ref = 14.864*scale^2; % m^2, aerodynamic reference area based 
            %on AMOOS reresentative packaging (18.3m long) original idea of X-37B is 3m high 
        in.v.aero.nose_radius = 3.7253*scale; % m, nose radius: AMOOS packaging (for 18.3 long)
        in.v.aero.cd = 3.874; % nd, cD at cL_max
        in.v.aero.cl = 2.003; % nd, cL_max
        in.v.aero.aoa = 45 * pi/180;    %trim angle of attack at cL_max

    %ESA IXV constant cL and cD
%         in.v.aero.mode = uint8(1); % nd, aerodynamics mode
%         in.v.aero.area_ref = 7.26; % m^2, aerodynamic reference area 
%         % in.v.aero.nose_radius = ; % m, nose radius:
%         in.v.aero.cd = 1.515; % nd, cD at cL_max
%         in.v.aero.cl = 1.0605; % nd, cL_max from published L/D value 
% %         in.v.aero.cl = 0.06;    %cL_max from data
%         in.v.aero.aoa = -25*pi/180;    %trim angle of attack at cL_max
        
% Decelerator
    %MSL DGB parachute - for 6m long sub
        in.v.decel.mode = uint8(2);       %mach dependent
        in.v.decel.d = 17;                %17m diameter DGB parachute
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

% Guidance
% %Bank
in.v.gnc.g.p.rate = 0.5; % Hz
in.v.gnc.g.p.bank_mode = uint8(2); % nd, constant bank rate
in.v.gnc.g.p.const_bank_rate = 4*pi/60; % rad/s, 2 RPM
in.v.gnc.g.p.const_bank = in.v.gnc.g.p.const_bank_rate; 
in.v.gnc.g.p.tgt_lat = tgt_lat; % rad, for range calculations
in.v.gnc.g.p.tgt_lon = tgt_lon; % rad, for range calculations

% Control
in.v.gnc.c.p.bank_mode = uint8(3);  %accel and rate-limit control
in.v.gnc.c.p.rate = 20; % Hz
in.v.gnc.c.p.decel_mode = uint8(1);
in.v.gnc.c.s.aoa = in.v.aero.aoa;
% in.v.gnc.c.p.aoa_mode = uint8(1);
in.v.gnc.c.s.bank_rate = 4*pi/60;

in.v.gnc.g.p.aoa_mode = uint8(1);
in.v.gnc.g.p.const_aoa = in.v.aero.aoa; 
in.s.traj.aoa_ini = in.v.aero.aoa;

%% Run Nominal Trajectory
mach_chute = 1.5;
[nom,nom2,mach_ind,ind2,vel,alt, ... 
    time,dr,cr,heat_rate,q,decel,lat,lon,mass,drag,az,gamma] = Titan_traj_slender(in,mach_chute,in.v.decel.d);

peak_heat = max(heat_rate)/(100^2);         %W/cm2   
heat_load = trapz(time,heat_rate)/(100^2);  %J/cm2
peak_g = max(decel);                        %Earth G's
peak_q = max(q);                            %Pa
descent_time = time(end)/3600;              %total descent time, hrs


%% Run Disperse Trajectories
in.s.traj.run_disp = 1;
if in.s.traj.run_disp == 1;
    n = 2000; %number of trajectories to run
    
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
    disp_ind2 = nan(1,n);                   %index when vehicle reaches the ground
    
    %initailize input variables for random sampling
    disp_mach_chute = nan(1,n);
    
    %other input variables
    tic;
    for i = 1:n
        i
        %apply randomly sampled uncertainties
        mc_val = randi(1000,1);     %generate random integer between 1 and 1000
        in_disp = apply_dispersions(in);
        in_disp.p.atm.table = atm_tables.fixed_atm_mc_tables.mc(mc_val).table;     %table of atmospheric data of
        %Titan-GRAM monte carlo simulation of random integer mc_val
        
        [in_disp.s.traj.r_pci_ini, in_disp.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
            in_disp.s.traj.lat, in_disp.s.traj.lon, in_disp.s.traj.alt, ...
            in_disp.s.traj.gamma_pp, in_disp.s.traj.az, in_disp.s.traj.vel_pp_mag, ...
            in_disp.p.omega, in_disp.s.traj.t_ini, 0, in_disp.p.r_e, in_disp.p.r_p);
        
        disp_mach_chute(i) = mach_chute + normrnd(0,0.01);    %mach number at chute deployment +- 0.01 1 sigma
        
        [~,~,disp_mach_ind(i),disp_ind2(i),disp_vel(:,i),disp_alt(:,i), ...
            disp_time(:,i),disp_dr(:,i),disp_cr(:,i),disp_heat_rate(:,i), disp_q(:,i),disp_decel(:,i), ...
            disp_lat(:,i),disp_lon(:,i),disp_mass(:,i),disp_drag(:,i),disp_az(:,i),disp_gamma(:,i)] = ...
            Titan_traj_slender(in_disp,disp_mach_chute(i),in.v.decel.d);
        
        %find end of trajectory
        sd_ind = find(isnan(disp_time(:,i)) == 1);
        sd_ind = sd_ind(1) - 1;
        
        v_sd(i) = disp_vel(sd_ind,i)';              %splashdown velocity
        disp_peak_heat(i) = max(disp_heat_rate(:,i))/(100^2)';        %peak heating rate in each case, W/cm2
        disp_heat_load(i) = trapz(disp_time(1:sd_ind,i),disp_heat_rate(1:sd_ind,i))/(100^2); %heat load in each case, J/cm2
        peak_g(i) = max(disp_decel(:,i))';           %peak deceleration, earth g's
        peak_q(i) = max(disp_q(:,i))';                   %peak dynamic pressure, N/m2
        descent_time(i) = disp_time(sd_ind,i)/3600';             %descent time each case, hrs
        dr_sd(i) = disp_dr(sd_ind,i)/1000';                  %downrange at splashdown. km
        cr_sd(i) = disp_cr(sd_ind,i)/1000';                  %crossrange at splashdown, km
        lat_sd(i) = disp_lat(sd_ind,i).*180/pi';                %latitude at splashdown, degrees
        lon_sd(i) = disp_lon(sd_ind,i).*180/pi';                %longitude at splashdown, degrees
        
        clear in_disp %clear data
    end
    toc;
    
    max_sd = max(v_sd);       %maximum splashdown speed
    min_sd = min(v_sd);       %mimimum splashdown speed
    avg_sd = mean(v_sd);      %average splashdown speed
    std_sd = std(v_sd);       %standard deviation splashdown speed
    med_sd = median(v_sd);
    
    %downrange stats
    max_dr = max(dr_sd);      %maximum downrange(km) from entry position
    min_dr = min(dr_sd);       %minimum downrange(km) from entry position
    avg_dr = mean(dr_sd);      %mean downrange(km) from entry position
    med_dr = median(dr_sd);
    dr_range = max_dr - min_dr;     %downrange landing ellipse
    
    %crossrange stats
    max_cr = max(cr_sd);      %maximum crossrange(km) from entry position
    min_cr = min(cr_sd);       %minimum crossrange(km) from entry position
    avg_cr = mean(cr_sd);      %mean crossrange(km) from entry position
    med_cr = median(cr_sd);
    cr_range = max_cr - min_cr;     %crossrange landing ellipse
    
    
    figure(1); hold on
    plot(lon_sd.*180/pi,lat_sd.*180/pi,'o','LineWidth',2,'MarkerSize',2)
    
    sea_map = imread(['.' filesep 'docs' filesep 'Ligeia_map.png']);
    imshow(sea_map);
    hold on
    axis off
   
%     axis([60 160 72 84])
%     imshow(sea_map);
%     truesize(sea_map)
    
%     set(gcf,'{pixels}','points','position',[0,0,1240,904])
    
    
    % figure(2); hold on;
    % plot(disp_vel,disp_alt);
    % plot(nom.traj.vel_pp_mag,nom.traj.alt,'LineWidth',1,'Color','k');
    % xlabel('Velocity (m/s)'); ylabel('Altitude (m)');
    % hold off
    %
    % figure(3); hold on
    % plot(disp_lat.* 180/pi,disp_lon.*180/pi);
    % plot(tgt_lat.*180/pi,tgt_lon.*180/pi,'*','LineWidth',3,'Color','k');
    % xlabel('Latitude  (deg)'); ylabel('Longitude (deg)');
    % set(gca,'FontSize',16)
    %
    % figure(4); hold on;
    % plot(disp_cr./1000, disp_alt);
    % plot(nom.traj.crossrange./1000,nom.traj.alt,'LineWidth',1,'Color','k');
    % xlabel('Crossrange (km)'); ylabel('Altitude (m)');
    % set(gca,'FontSize',16)
    %
    % figure(5); hold on;
    % plot(disp_dr./1000,disp_alt);
    % plot(nom.traj.downrange./1000,nom.traj.alt,'LineWidth',1,'Color','k');
    % xlabel('Downrange (km)'); ylabel('Altitude (m)');
    % set(gca,'FontSize',16)
    %
    % figure(6); hold on;
    % plot(disp_dr./1000,disp_cr./1000);
    % xlabel('Downrange (km)'); ylabel('Crossrange (km)');
    % set(gca,'FontSize',16)
end