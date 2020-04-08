% Titan Submarine
%   Trajectory Analysis for Titan Turtle submarine
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
%
%
% Major Revision History:
%   *Created SEP 2016 - E. ROELKE
%#codegen
clear; 
clc;
% clf;
format longg;
%% Define input structure
    in = default_in;

%% Planetary data
    in.p = get_planetary_data(6,3); % Load settings for Titan (6), model with winds (3)
    
    winds = 'sin';
    if strcmp(winds,'sin') == 1
        atm_tables = load([ '.' filesep 'data' filesep 'atm_titan_gram_10000mc_sin_winds.mat']);
        in.p.atm.table = atm_tables.atm_mc_sin_tables.nom_table;
    else
        atm_tables = load([ '.' filesep 'data' filesep 'atm_titan_gram_10000mc_fixed_winds.mat']);
        in.p.atm.table = atm_tables.fixed_atm_mc_tables.nom_table; winds = 'regular';
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
    in.s.traj.t_max = 15000; % double, s, maximum allowable time
    in.s.data_rate = 1;

% Initial state (Eastbound equatorial)
    in.s.traj.alt = 1270e3; % atmospheric interface altitude, m
    in.s.traj.vel_pp_mag = 4000; % initial velocity, m/s
    
    in.s.traj.lat = 81 * pi/180; % initial latitude (72.4deg north), rad
    in.s.traj.lon = -338*pi/180; % initial longitude (338deg west), rad
    in.s.traj.gamma_pp = -65.4*pi/180; % initial flight-path, rad (Labreton ref)
    in.s.traj.az = 88.89*pi/180; % initial azimuth, rad

% Convert scalars to inertial position and velocity vectors
    [in.s.traj.r_pci_ini, in.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
        in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
        in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
        in.p.omega, in.s.traj.t_ini, 0, in.p.r_e, in.p.r_p);

%% Vehicle data
% Mass properties
    in.v.mp.m_ini = 600; % kg, vehicle mass at entry

% Aerodynamics - mach dependent cl and cd for ballistic 70deg sphere cone
    in.v.aero.mode = uint8(2); % nd, table look up
    temp = load( ['.' filesep 'data' filesep 'aero_sphere_cone_70deg.mat'] );
    in.v.aero.table = temp.table;
    in.v.aero.nose_radius = 1.25;  %nose radius, meters
    in.v.aero.area_ref = pi*(2.7/2)^2;  %reference area, 2.7m diam aeroshell
    d = 2*(in.v.aero.area_ref/pi)^(1/2);        %diameter of aeroshell
    in.v.aero.chord_ref = 1.97;      %1.9m long aeroshell
    in.v.aero.aoa = 0; %angle of attack, radians

% Decelerator
    %DGB Pilot chute for sphere cone aeroshell
        in.v.decel.mode = uint8(1); %user-input
        in.v.decel.d = 3; % parachute diameter, m
        in.v.decel.K = 0; % no dispersions (nominal)
        in.v.decel.cd = 0.55;

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

    %deployment parameters
        mach_deploy = 1.5;
        v_deploy = 20;
        d_pf = 7;
        d_pc = 3;
%% Trajectory Analysis
    [nom,nom2,nom3,terminate,mach_ind,vel_ind2,ind3,vel,alt,time,dr,cr,heat_rate,q,decel,lat,lon,mass,drag,az,gamma] = ...
        TSSM_Traj(in,mach_deploy,v_deploy,d_pf,d_pc,in.s.traj.run_disp);
    
    %{
    %plots
    %split up trajectories
        figure(1); hold on
        set(gca,'FontSize',16)
        plot(nom.traj.vel_pp_mag(1:mach_ind),nom.traj.alt(1:mach_ind)./1000,'b')
        plot(nom2.traj.vel_pp_mag(1:vel_ind2),nom2.traj.alt(1:vel_ind2)./1000,'r')
        plot(nom3.traj.vel_pp_mag,nom3.traj.alt./1000,'g')
        xlabel('Planet-Relative velocity, m/s')
        ylabel('Altitude, km')
        legend('Hypersonic Descent','Pilot Chute Deployed','Parafoil Deployed','Location','best');
% 
        figure(2); hold on
        set(gca,'FontSize',16)
        plot(nom.traj.time(1:mach_ind)./3600,nom.traj.alt(1:mach_ind)./1000,'b')
        plot(nom2.traj.time(1:vel_ind2)./3600,nom2.traj.alt(1:vel_ind2)./1000,'r')
        plot(nom3.traj.time./3600,nom3.traj.alt./1000,'g')
        xlabel('Time, hrs')
        ylabel('Altitude, km')
        legend('Hypersonic Descent','Pilot Chute Deployed','Parafoil Deployed','Location','best');
%     
        %combined trajectories
        figure(3);
        set(gca,'FontSize',16)
        plot(vel,alt)
        xlabel('Planet-Relative velocity, m/s')
        ylabel('Altitude, km')
    %}

%find peak values
    peak_heat = max(heat_rate)/(100^2);        %peak heating rate, W/cm2
    peak_g = max(decel);           %peak deceleration, earth g's
    peak_q = max(q);                   %peak dynamic pressure, N/m2

%total trajectory values
    descent_time = time(end)/3600;     %descent time, hrs

    heat_load = trapz(time,heat_rate)/(100^2);  %J/cm2

    if alt(end) < 5
        v_splashdown = nom3.traj.vel_pp_mag(ind3);      %splashdown velocity
    else
        disp('vehicle did not reach the ground')
        nom3.traj.alt(ind3)
    end


    
%% Disperse Trajectories
% %{
in.s.traj.run_disp = 0;
if in.s.traj.run_disp == 1;
    %run monte carlo / disperse trajectories
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
        disp_vel_ind2 = nan(1,n);                   %index when parafoil deploys
        disp_ind3 = nan(1,n);                       %index when vehicle reaches the ground
    
        %initailize input variables for random sampling
        disp_mach_pc = nan(1,n);
        disp_v_pf = nan(1,n);
        disp_terminate = nan(1,n);

%other input variables
tic;
for i = 1:n
    i
    
    %apply randomly sampled uncertainties
    mc_val = randi(1000,1);     %generate random integer between 1 and 1000
    in_disp = apply_dispersions(in); 
    if strcmp(winds,'sin') == 1
        in_disp.p.atm.table = atm_tables.atm_mc_sin_tables.mc(mc_val).table;
    else
        in_disp.p.atm.table = atm_tables.fixed_atm_mc_tables.mc(mc_val).table;
    end
        
        
    [in_disp.s.traj.r_pci_ini, in_disp.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
        in_disp.s.traj.lat, in_disp.s.traj.lon, in_disp.s.traj.alt, ...
        in_disp.s.traj.gamma_pp, in_disp.s.traj.az, in_disp.s.traj.vel_pp_mag, ...
        in_disp.p.omega, in_disp.s.traj.t_ini, 0, in_disp.p.r_e, in_disp.p.r_p);
    
    disp_mach_pc(i) = mach_deploy + normrnd(0,0.01);    %mach number at pilot chute deployment +- 0.01 1 sigma
    disp_v_pf(i) = v_deploy + normrnd(0,1);             %velocity at parafoil deployment +- 1 m/s 1 sigma
   
    [~,~,~,disp_terminate(i),disp_mach_ind(i),disp_vel_ind2(i),disp_ind3(i),disp_vel(:,i),disp_alt(:,i), ... 
    disp_time(:,i),disp_dr(:,i),disp_cr(:,i),disp_heat_rate(:,i), disp_q(:,i),disp_decel(:,i), ... 
    disp_lat(:,i),disp_lon(:,i),disp_mass(:,i),disp_drag(:,i),disp_az(:,i),disp_gamma(:,i)] = ... 
    TSSM_Traj(in_disp,disp_mach_pc(i),disp_v_pf(i),d_pf,d_pc,in.s.traj.run_disp);
    
    %find end of trajectory
    sd_ind = find(isnan(disp_time(:,i)) == 1);
    sd_ind = sd_ind(1) - 1;

    v_sd(i) = disp_vel(sd_ind,i)';              %splashdown velocity
    peak_heat(i) = max(disp_heat_rate(:,i))/(100^2)';        %peak heating rate, W/cm2
    heat_load(i) = trapz(disp_time(1:sd_ind,i),disp_heat_rate(1:sd_ind,i)); %heat load in each case
    peak_g(i) = max(disp_decel(:,i))';           %peak deceleration, earth g's
    peak_q(i) = max(disp_q(:,i))';                   %peak dynamic pressure, N/m2
    descent_time(i) = disp_time(sd_ind,i)/3600';             %descent time each case, hrs
    dr_sd(i) = disp_dr(sd_ind,i)/1000';                  %downrange at splashdown. km
    cr_sd(i) = disp_cr(sd_ind,i)/1000';                  %crossrange at splashdown, km
    lat_sd(i) = disp_lat(sd_ind,i)';                %latitude at splashdown, rad
    lon_sd(i) = disp_lon(sd_ind,i)';                %longitude at splashdown, rad
    if disp_alt(sd_ind) > 1
        disp('vehicle did not reach the ground')
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
    avg_dr = mean(dr_sd);      %mean downrange(km) from entry position
    med_dr = median(dr_sd);     %median downrange(km) from entry
    dr_range = max_dr - min_dr;     %downrange landing ellipse

%crossrange stats
    max_cr = max(cr_sd);      %maximum crossrange(km) from entry position
    min_cr = min(cr_sd);       %minimum crossrange(km) from entry position
    avg_cr = mean(cr_sd);      %mean crossrange(km) from entry position
    med_cr = median(cr_sd);     %median crossrange(km) from entry
    cr_range = max_cr - min_cr;     %crossrange landing ellipse
    
%latitude, longitude
    lat_sd = lat_sd.*180/pi;            %degrees N
    lon_sd = 360 - lon_sd.*180/pi;      %degrees W
    
    
figure(1); hold on
plot(lon_sd,lat_sd,'o','MarkerSize',2,'LineWidth',2)
xlabel('Longitude W'); ylabel('Latitude N')
   

figure(2); hold on;
plot(disp_vel,disp_alt);
xlabel('Planet-Relative velocity, m/s')
ylabel('Altitude, km')
hold off
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



%% Parachute Size Trade Study
trade_study = 0;
if trade_study == 1
    clear; clc;
    D_steerable = 3:0.5:15; %main parachute diameters, m
    Mach_main = nom.traj.mach(mach_ind);    %mach number at pilot chute deployment
    D_pilot = in2.v.decel.d;            %pilot chute diameter
    V_steerable = v_pf;       %m/s
    for i = 1:length(D_steerable)
        [descent_time(i), heat_load(i), v_impact(i)] = ...
            titan_parachute_tradestudy(D_steerable(i),Mach_main,D_pilot,V_steerable);
    end
    
    figure(1)
    plot(descent_time,D_steerable)
    set(gca,'FontSize',16)
    xlabel('Descent Time (hrs)')
    ylabel('Main Parachute Diameter, m')

    figure(2)
    plot(v_impact,D_steerable)
    xlabel('Impact velocity, m/s')
    ylabel('Main Parachute Diameter, m')
end
