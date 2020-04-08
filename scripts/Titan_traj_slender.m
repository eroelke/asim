% 3 phase trajectory for titan EDL
%   1. ballistic
%   2. pilot parachute
%   3. main parachute
%   
% Written by:
%   E. Roelke, Spring 2016
% 
function [nom,nom2,mach_ind,ind2,vel,alt, ... 
    time,dr,cr,heat_rate,q,decel,lat,lon,mass,drag,az,gamma] = Titan_traj_slender(in,mach_chute,d_chute)
%#codegen

in.v.decel.d = d_chute;

nom = main1_mex(in);

%find index and important values (nominal) when vehicle reaches mach_chute
mach_ind = find(nom.traj.mach <= mach_chute);
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

%DGB parachute deployment
    in2 = in;          
    in2.s.traj.t_max = 20000;
    in2.s.traj.t_ini = nom.traj.time(mach_ind);
    in2.s.traj.alt = nom.traj.alt(mach_ind);     		%parachute opens
    in2.s.traj.vel_pp_mag = nom.traj.vel_pp_mag(mach_ind);    	%planet-relative velocity magnitude
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

    in2.v.gnc.c.p.decel_mode = uint8(1);
    in2.v.gnc.c.s.deploy_decel = logical(true);
    in2.v.gnc.g.p.prop_mode = uint8(10);        %split up trajectories, declerator already deployed

    nom2 = main1_mex(in2);
    ind2 = find(isnan(nom2.traj.time) == 1);
    ind2 = ind2(1) - 1;
    
    
    
%% Combine Data
    if in.s.traj.run_disp == 0
        vel = vertcat(nom.traj.vel_pp_mag(1:mach_ind),nom2.traj.vel_pp_mag(2:ind2));
        alt = vertcat(nom.traj.alt(1:mach_ind),nom2.traj.alt(2:ind2));
        time = vertcat(nom.traj.time(1:mach_ind),nom2.traj.time(2:ind2));
        dr = vertcat(nom.traj.downrange(1:mach_ind),nom2.traj.downrange(2:ind2));
        cr = vertcat(nom.traj.crossrange(1:mach_ind),nom2.traj.crossrange(2:ind2));
        heat_rate = vertcat(nom.traj.heat_rate(1:mach_ind),nom2.traj.heat_rate(2:ind2));
        q = vertcat(nom.traj.q(1:mach_ind),nom2.traj.q(2:ind2));
        decel = vertcat(nom.traj.g_loading(1:mach_ind),nom2.traj.g_loading(2:ind2));
        lat = vertcat(nom.traj.lat(1:mach_ind),nom2.traj.lat(2:ind2));
        lon = vertcat(nom.traj.lon(1:mach_ind),nom2.traj.lon(2:ind2));
        mass = vertcat(nom.traj.mass(1:mach_ind),nom2.traj.mass(2:ind2));
        drag = vertcat(nom.traj.drag_ii_mag(1:mach_ind),nom2.traj.drag_ii_mag(2:ind2));
        az = vertcat(nom.traj.azi_pp(1:mach_ind),nom2.traj.azi_pp(2:ind2));
        gamma = vertcat(nom.traj.gamma_pp(1:mach_ind),nom2.traj.gamma_pp(2:ind2));
    elseif in.s.traj.run_disp == 1      %need to include NaN to keep column length the same
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
    else
        error('in.s.traj.run_disp must be 0 or 1');
    end
    
    
%check to see if vehicle reaches ground
if nom2.traj.alt(ind2) < 5
    v_splashdown = vel(ind2);      %splashdown velocity
else
    disp('vehicle did not reach the ground')
    alt(ind2)
end

end
