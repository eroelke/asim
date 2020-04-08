% 3 phase trajectory for titan EDL
%   1. ballistic
%   2. pilot parachute
%   3. main parachute
%   
% Written by:
%   E. Roelke, Spring 2016
function [nom,nom2,nom3,terminate,mach_ind,vel_ind2,ind3,vel,alt, ...
    time,dr,cr,heat_rate,q,decel,lat,lon,mass,drag,az,gamma] = TitanTraj_TwoChutes(in,mach_pc,v_pf,d_pf,d_pc)
%#codegen

terminate = 0;
%pilot chute
in.v.decel.d = d_pc;
%equations of motion integration
nom = main1_mex(in);

%find index and important values (nominal) when vehicle reaches mach 1.5
mach_ind = find(nom.traj.mach <= mach_pc);
mach_ind = mach_ind(1);

%pilot parachute deployment
in2 = in;
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

%find index and important values (nominal) when vehicle reaches
%v_pf
temp = find(nom2.traj.vel_pp_mag <= v_pf);
if isempty(temp) == 1
    terminate = 1;
    vel_ind2 = 0; nom3 = 0; ind3 = 0;
    temp = find(isnan(nom2.traj.vel_pp_mag) == 1) - 1; temp = temp(1);
    %             v_splashdown = nom2.traj.vel_pp_mag(temp);
    %             keyboard
else
    vel_ind2 = temp(1);
end

if terminate == 0
    %parafoil deployment
    in3 = in;
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
    %if dispersion analysis, apply dispersion to drag coefficient
    if in.s.traj.run_disp == 1
        in3.v.decel.cd = in3.v.decel.cd + normrnd(0,0.03);
    end
    in3.v.decel.d = d_pf;                             %m, parafoil diameter
    in3.v.mp.m_ini = in2.v.mp.m_ini - 78;           %eject front/backshell, pilot chute, TPS
    in3.v.gnc.c.s.deploy_decel = logical(true);
    in3.v.gnc.g.p.prop_mode = uint8(10);             %split up trajectories, declerator already deployed
    
    nom3 = main1_mex(in3);
    ind3 = find(isnan(nom3.traj.time) == 1);
    ind3 = ind3(1) - 1;
    
    
    %post-processing
    %combine data
    if in.s.traj.run_disp == 0
        vel = vertcat(nom.traj.vel_pp_mag(1:mach_ind),nom2.traj.vel_pp_mag(2:vel_ind2),nom3.traj.vel_pp_mag(2:ind3));
        alt = vertcat(nom.traj.alt(1:mach_ind),nom2.traj.alt(2:vel_ind2),nom3.traj.alt(2:ind3));
        time = vertcat(nom.traj.time(1:mach_ind),nom2.traj.time(2:vel_ind2),nom3.traj.time(2:ind3));
        dr = vertcat(nom.traj.downrange(1:mach_ind),nom2.traj.downrange(2:vel_ind2),nom3.traj.downrange(2:ind3));
        cr = vertcat(nom.traj.crossrange(1:mach_ind),nom2.traj.crossrange(2:vel_ind2),nom3.traj.vel_pp_mag(2:ind3));
        heat_rate = vertcat(nom.traj.heat_rate(1:mach_ind),nom2.traj.heat_rate(2:vel_ind2),nom3.traj.heat_rate(2:ind3));
        q = vertcat(nom.traj.q(1:mach_ind),nom2.traj.q(2:vel_ind2),nom3.traj.q(2:ind3));
        decel = vertcat(nom.traj.g_loading(1:mach_ind),nom2.traj.g_loading(2:vel_ind2),nom3.traj.g_loading(2:ind3));
        lat = vertcat(nom.traj.lat(1:mach_ind),nom2.traj.lat(2:vel_ind2),nom3.traj.lat(2:ind3));
        lon = vertcat(nom.traj.lon(1:mach_ind),nom2.traj.lon(2:vel_ind2),nom3.traj.lon(2:ind3));
        mass = vertcat(nom.traj.mass(1:mach_ind),nom2.traj.mass(2:vel_ind2),nom3.traj.mass(2:ind3));
        drag = vertcat(nom.traj.drag_ii_mag(1:mach_ind),nom2.traj.drag_ii_mag(2:vel_ind2),nom3.traj.drag_ii_mag(2:ind3));
        az = vertcat(nom.traj.azi_pp(1:mach_ind),nom2.traj.azi_pp(2:vel_ind2),nom3.traj.azi_pp(2:ind3));
        gamma = vertcat(nom.traj.gamma_pp(1:mach_ind),nom2.traj.gamma_pp(2:vel_ind2),nom3.traj.gamma_pp(2:ind3));
    elseif in.s.traj.run_disp == 1
        vel = vertcat(nom.traj.vel_pp_mag(1:mach_ind),nom2.traj.vel_pp_mag(2:vel_ind2),nom3.traj.vel_pp_mag(2:end));
        alt = vertcat(nom.traj.alt(1:mach_ind),nom2.traj.alt(2:vel_ind2),nom3.traj.alt(2:end));
        time = vertcat(nom.traj.time(1:mach_ind),nom2.traj.time(2:vel_ind2),nom3.traj.time(2:end));
        dr = vertcat(nom.traj.downrange(1:mach_ind),nom2.traj.downrange(2:vel_ind2),nom3.traj.downrange(2:end));
        cr = vertcat(nom.traj.crossrange(1:mach_ind),nom2.traj.crossrange(2:vel_ind2),nom3.traj.vel_pp_mag(2:end));
        heat_rate = vertcat(nom.traj.heat_rate(1:mach_ind),nom2.traj.heat_rate(2:vel_ind2),nom3.traj.heat_rate(2:end));
        q = vertcat(nom.traj.q(1:mach_ind),nom2.traj.q(2:vel_ind2),nom3.traj.q(2:end));
        decel = vertcat(nom.traj.g_loading(1:mach_ind),nom2.traj.g_loading(2:vel_ind2),nom3.traj.g_loading(2:end));
        lat = vertcat(nom.traj.lat(1:mach_ind),nom2.traj.lat(2:vel_ind2),nom3.traj.lat(2:end));
        lon = vertcat(nom.traj.lon(1:mach_ind),nom2.traj.lon(2:vel_ind2),nom3.traj.lon(2:end));
        mass = vertcat(nom.traj.mass(1:mach_ind),nom2.traj.mass(2:vel_ind2),nom3.traj.mass(2:end));
        drag = vertcat(nom.traj.drag_ii_mag(1:mach_ind),nom2.traj.drag_ii_mag(2:vel_ind2),nom3.traj.drag_ii_mag(2:end));
        az = vertcat(nom.traj.azi_pp(1:mach_ind),nom2.traj.azi_pp(2:vel_ind2),nom3.traj.azi_pp(2:end));
        gamma = vertcat(nom.traj.gamma_pp(1:mach_ind),nom2.traj.gamma_pp(2:vel_ind2),nom3.traj.gamma_pp(2:end));
    end
elseif terminate == 1
    if in.s.traj.run_disp == 0
        vel = vertcat(nom.traj.vel_pp_mag(1:mach_ind),nom2.traj.vel_pp_mag(2:temp));
        alt = vertcat(nom.traj.alt(1:mach_ind),nom2.traj.alt(2:temp));
        time = vertcat(nom.traj.time(1:mach_ind),nom2.traj.time(2:temp));
        dr = vertcat(nom.traj.downrange(1:mach_ind),nom2.traj.downrange(2:temp));
        cr = vertcat(nom.traj.crossrange(1:mach_ind),nom2.traj.crossrange(2:temp));
        heat_rate = vertcat(nom.traj.heat_rate(1:mach_ind),nom2.traj.heat_rate(2:temp));
        q = vertcat(nom.traj.q(1:mach_ind),nom2.traj.q(2:temp));
        decel = vertcat(nom.traj.g_loading(1:mach_ind),nom2.traj.g_loading(2:temp));
        lat = vertcat(nom.traj.lat(1:mach_ind),nom2.traj.lat(2:temp));
        lon = vertcat(nom.traj.lon(1:mach_ind),nom2.traj.lon(2:temp));
        mass = vertcat(nom.traj.mass(1:mach_ind),nom2.traj.mass(2:temp));
        drag = vertcat(nom.traj.drag_ii_mag(1:mach_ind),nom2.traj.drag_ii_mag(2:temp));
        az = vertcat(nom.traj.azi_pp(1:mach_ind),nom2.traj.azi_pp(2:temp));
        gamma = vertcat(nom.traj.gamma_pp(1:mach_ind),nom2.traj.gamma_pp(2:temp));
    elseif in.s.traj.run_disp == 1
        vel = vertcat(nom.traj.vel_pp_mag(1:mach_ind),nom2.traj.vel_pp_mag(2:end));
        alt = vertcat(nom.traj.alt(1:mach_ind),nom2.traj.alt(2:end));
        time = vertcat(nom.traj.time(1:mach_ind),nom2.traj.time(2:end));
        dr = vertcat(nom.traj.downrange(1:mach_ind),nom2.traj.downrange(2:end));
        cr = vertcat(nom.traj.crossrange(1:mach_ind),nom2.traj.crossrange(2:end));
        heat_rate = vertcat(nom.traj.heat_rate(1:mach_ind),nom2.traj.heat_rate(2:end));
        q = vertcat(nom.traj.q(1:mach_ind),nom2.traj.q(2:end));
        decel = vertcat(nom.traj.g_loading(1:mach_ind),nom2.traj.g_loading(2:end));
        lat = vertcat(nom.traj.lat(1:mach_ind),nom2.traj.lat(2:end));
        lon = vertcat(nom.traj.lon(1:mach_ind),nom2.traj.lon(2:end));
        mass = vertcat(nom.traj.mass(1:mach_ind),nom2.traj.mass(2:end));
        drag = vertcat(nom.traj.drag_ii_mag(1:mach_ind),nom2.traj.drag_ii_mag(2:end));
        az = vertcat(nom.traj.azi_pp(1:mach_ind),nom2.traj.azi_pp(2:end));
        gamma = vertcat(nom.traj.gamma_pp(1:mach_ind),nom2.traj.gamma_pp(2:end));
    end
end

% %find peak values
%     peak_heat = max(heat_rate)/(100^2);        %peak heating rate, W/cm2
%     peak_g = max(decel);           %peak deceleration, earth g's
%     peak_q = max(q);                   %peak dynamic pressure, N/m2
end