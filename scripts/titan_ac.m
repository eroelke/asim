%% titan_ac
clear;clc





% TCW plots
efpa0 = -30;
v0s = 5:1:12;
betas = [50 100 150 200 500];

[x0,aero,gnc,sim,mc] = base_jsr();
sim.planet = 'titan';

% Ni = length(efpas);
Nj = length(v0s);
Nk = length(betas);
fpa = nan(Nj,Nk);

cd = 1.2;
rc = 4.5/2;
m = nan(Nk,1);

fpa_guess = [-5 -50];
ha_tgt = 2000;
tol_fpa = 0.2;
tol_ha = 10;

in = dme_titan_in;
in.v.aero.area_ref = pi * rc^2;
in.v.aero.cd = cd;
in.v.aero.nose_radius = 0.5;
in.v.aero.mode = uint8(1);
in.v.gnc.g.p.rate = 0;

in.s.traj.alt = 1270e3; % atmospheric interface altitude, m
in.s.traj.lat = 73.4 * pi/180; % initial latitude (72.4deg north), rad
in.s.traj.lon = -338*pi/180; % initial longitude (338deg west), rad
% in.s.traj.vel_pp_mag = 6000; % initial velocity, m/s
% in.s.traj.gamma_pp = -65.4*pi/180; % initial flight-path, rad (Labreton ref)
in.s.traj.gamma_pp = efpa0 * pi/180;
in.s.traj.az = 88.89*pi/180; % initial azimuth, rad
in.s.traj.t_max = 8000;
in.s.traj.rate = 100;

for k = 1:Nk
    for j = 1:Nj
            
        m(k) = betas(k) * pi * cd * rc^2;
        in.v.mp.m_ini = m(k);
        in.s.traj.vel_pp_mag = v0s(j)*1000;

        % Convert scalars to inertial position and velocity vectors
        [in.s.traj.r_pci_ini, in.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
            in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
            in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
            in.p.omega, in.s.traj.t_ini, 0, in.p.r_e, in.p.r_p);

            fpa(j,k) = find_fpa_given_ha( ...
                fpa_guess, ... % 1x2 vector of initial guess on bounds, MUST encompass solution (deg)
                ha_tgt, ... % target apoapse altitude (km)
                in, ... % asim input file, implicitly contains plantary parameters and ICs
                tol_fpa, ... % tolerance on final output FPA (deg)
                tol_ha); % tolerance on final target apoapse altitude (km)


%         [fpa(j,k),~,haf_err,~] = fpa_given_haf( ...
%             efpa0, ... %  initial fpa guess, (deg)
%             ha_tgt, ... % target apoapse altitude, km)
%             in, ... % asim input file, contains plantary parameters, vehicle, and ICs
%             tol_fpa, ... % tolerance on flight path angle to search for (deg)
%             tol_ha); % tolerance on final target apoapse altitude km)
        
    end
end