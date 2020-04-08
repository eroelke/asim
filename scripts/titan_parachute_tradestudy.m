%Titan Parachute Trade Study
    %Perform a trade study on parachute diameter(s) and opening velocities as
    %a function of descent time
    
    %Inputs:
        %D_main: diameter of main parachute
        %Mach_main: mach number at which main parachute opens
        %D_steerable: diameter of guided parachute
        %V_steerable: velocity at which to open guided chute
    
    %Outputs:
        %descent_time: total trajectory descent time
        %
% 
% Written by:
%   E. Roelke, Spring 2016
% 
function [descent_time, heat_load, v_impact, alt_pc, alt_pf, t_pc, t_pf] = ... 
    titan_parachute_tradestudy(D_steerable,Mach_main,D_pilot,V_steerable)
%#codegen
   in = default_in;

%% Planetary data
    %TitanGRAM may be innacurate with wind calcs, especially at low
    %altitudes
    in.p = get_planetary_data(6,3); % Load settings for Titan (6), model with winds (3)
    temp = load([ '.' filesep 'data' filesep 'atm_titan_gram_nominal.mat']);    %temp table of atmospheric data
    in.p.atm.table = temp.table;


%% Simulation data
% Termination conditions
    in.s.term.or.type(1) = int8(0); % int8, nd, crossing/equality condition
    in.s.term.or.var(1) = uint8(1); % uint8, nd, termination mode (altitude)
    in.s.term.or.value(1) = 0; % double, m, terminate at 0 km altitude (sea level)
    in.s.term.or.direction(1) = int8(-1); % int8, nd, negative crossing

% Integration parameters
    in.s.traj.rate = 10; % double, Hz, integration rate
    in.s.traj.t_ini = 0; % double, s, initial time
    in.s.traj.t_max = 8000; % double, s, maximum allowable time
    in.s.data_rate = 1;

% Initial state (Eastbound equatorial)
    in.s.traj.alt = 1270e3; % atmospheric interface altitude, m
    in.s.traj.vel_pp_mag = 3000; % initial velocity, m/s
    
    in.s.traj.lat = 73.4 * pi/180; % initial latitude (72.4deg north), rad
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
        in.v.decel.d = D_pilot; % pilot chute diameter, m
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
    in.v.gnc.g.p.tgt_lat = 78*pi/180; % rad (70deg N - Kraken Mare center)
    in.v.gnc.g.p.tgt_lon = 45*pi/180;        %315 W = 45E

% Control
    in.v.gnc.c.p.bank_mode = uint8(3);  %accel and rate-limit control
    in.v.gnc.c.p.rate = 20; % Hz
    in.v.gnc.c.p.decel_mode = uint8(1);
    in.v.gnc.c.s.aoa = in.v.aero.aoa;
    in.v.gnc.c.s.bank_rate = 4*pi/60;

    in.v.gnc.g.p.aoa_mode = uint8(1);
    in.v.gnc.g.p.const_aoa = in.v.aero.aoa; 
    in.s.traj.aoa_ini = in.v.aero.aoa;


%% Trajectory Analysis
%equations of motion integration
    nom = main1_mex(in);
    ind = find(isnan(nom.traj.time) == 1);
    ind = ind(1) - 1;

    %find index and important values (nominal) when vehicle reaches mach 1.5
    mach_ind = find(nom.traj.mach <= Mach_main);
    mach_ind = mach_ind(1);
    v_pc = nom.traj.vel_pp_mag(mach_ind);    %vehicle velocity at pilot chute deployment
    alt_pc = nom.traj.alt(mach_ind);         %altitude at pilot chute deployment
    q_pc = nom.traj.q(mach_ind);            %dynamic pressure at pilot chute deployment
    t_pc = nom.traj.time(mach_ind)/60;         %time pilot chute deployed, min
    lat_pc = nom.traj.lat(mach_ind);        
    lon_pc = nom.traj.lon(mach_ind);
    cr_pc = nom.traj.crossrange(mach_ind);
    dr_pc = nom.traj.downrange(mach_ind);
    az_pc = nom.traj.azi_pp(mach_ind);


%pilot parachute deployment
    in2 = in;          
    in2.s.traj.t_max = 15000;                   %increase max time due to long descent time
    in2.s.traj.t_ini = nom.traj.time(mach_ind);
    in2.s.traj.alt = nom.traj.alt(mach_ind);     %parachute opens at 160km nominal
    in2.s.traj.vel_pp_mag = nom.traj.vel_pp_mag(mach_ind);    %planet-relative velocity magnitude
    in2.s.traj.az = nom.traj.azi_pp(mach_ind);
    in2.s.traj.cr_ini = nom.traj.crossrange(mach_ind);
    in2.s.traj.dr_ini = nom.traj.downrange(mach_ind);
    in2.s.traj.gamma_pp = nom.traj.gamma_pp(mach_ind);
    in2.s.traj.bank_ini = nom.traj.bank(mach_ind);
    in2.s.traj.aoa_ini = nom.traj.aoa(mach_ind);
    [in2.s.traj.r_pci_ini, in2.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
        in2.s.traj.lat, in2.s.traj.lon, in2.s.traj.alt, ...
        in2.s.traj.gamma_pp, in2.s.traj.az, in2.s.traj.vel_pp_mag, ...
        in2.p.omega, in2.s.traj.t_ini, 0, in2.p.r_e, in2.p.r_p);

    %Bank Control - turn off banking
    in2.v.gnc.g.p.const_bank_rate = 0; % rad/s
    
    %Turn Deploy_Decel on
    in2.v.gnc.c.s.deploy_decel = logical(true);
    in2.v.gnc.g.p.prop_mode = uint8(10);        %split up trajectories, declerator already deployed

    %Trajectory Analysis
        nom2 = main1_mex(in2);
        ind2 = find(isnan(nom2.traj.time) == 1);
        ind2 = ind2(1) - 1;

        %find index and important values (nominal) when vehicle reaches 40 m/s
        vel_ind2 = find(nom2.traj.vel_pp_mag <= V_steerable);   
        vel_ind2 = vel_ind2(1);
        v_pf = nom2.traj.vel_pp_mag(vel_ind2);    %velocity at parafoil deployment
        alt_pf = nom2.traj.alt(vel_ind2);         %altitude at parafoil deployment
        q_pf = nom2.traj.q(vel_ind2);            %dynamic pressure at parafoil deployment
        t_pf = nom2.traj.time(vel_ind2)/60;         %time at parafoil deployment, min
        lat_pf = nom2.traj.lat(vel_ind2);        
        lon_pf = nom2.traj.lon(vel_ind2);
        cr_pf = nom2.traj.crossrange(vel_ind2);
        dr_pf = nom2.traj.downrange(vel_ind2);
        az_pf = nom2.traj.azi_pp(vel_ind2);

%parafoil deployment
    in3 = in;           
    in3.s.traj.t_max = 20000;
    in3.s.traj.t_ini = nom2.traj.time(vel_ind2);
    in3.s.traj.alt = nom2.traj.alt(vel_ind2);     
    in3.s.traj.vel_pp_mag = nom2.traj.vel_pp_mag(vel_ind2);    %planet-relative velocity magnitude
    in3.s.traj.az = nom2.traj.azi_pp(vel_ind2);
    in3.s.traj.cr_ini = nom2.traj.crossrange(vel_ind2);
    in3.s.traj.dr_ini = nom2.traj.downrange(vel_ind2);
    in3.s.traj.gamma_pp = nom2.traj.gamma_pp(vel_ind2);
    in3.s.traj.bank_ini = nom2.traj.bank(vel_ind2);
    in3.s.traj.aoa_ini = nom2.traj.aoa(vel_ind2);
    [in3.s.traj.r_pci_ini, in3.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
        in3.s.traj.lat, in3.s.traj.lon, in3.s.traj.alt, ...
        in3.s.traj.gamma_pp, in3.s.traj.az, in3.s.traj.vel_pp_mag, ...
        in3.p.omega, in3.s.traj.t_ini, 0, in3.p.r_e, in3.p.r_p);

    in3.v.gnc.g.p.const_bank_rate = 0; % rad/s
        
    in3.v.gnc.c.p.decel_mode = uint8(1);
    in3.v.decel.cd = 0.55;
    in3.v.decel.d = D_steerable;                             %m, parafoil diameter
    in3.v.mp.m_ini = in2.v.mp.m_ini - 185;           %eject front/backshell, pilot chute, TPS
    in3.v.gnc.c.s.deploy_decel = logical(true);
    in3.v.gnc.g.p.prop_mode = uint8(10);             %split up trajectories, declerator already deployed

    nom3 = main1_mex(in3);
    ind3 = find(isnan(nom3.traj.time) == 1);
    ind3 = ind3(1) - 1;
    
%find peak values
    peak_heat = max(nom.traj.heat_rate)/(100^2);        %peak heating rate, W/cm2
    peak_g = max(nom.traj.g_loading);           %peak deceleration, earth g's
    peak_q = max(nom.traj.q);                   %peak dynamic pressure, N/m2

%total trajectory values
    descent_time = nom3.traj.time(ind3) / 3600;     %descent time, hrs

    heat_load = (1/100^2)*( trapz(nom.traj.time(1:mach_ind),nom.traj.heat_rate(1:mach_ind)) + ...
        trapz(nom2.traj.time(1:vel_ind2),nom2.traj.heat_rate(1:vel_ind2)) + ...
        trapz(nom3.traj.time(1:ind3),nom3.traj.heat_rate(1:ind3)) );            %J/cm2

    if nom3.traj.alt(ind3) < 5
        v_impact = nom3.traj.vel_pp_mag(ind3);      %splashdown velocity
    else
        disp('vehicle did not reach the ground')
        nom3.traj.alt(ind3)
        v_impact = nan;
    end

    
end

