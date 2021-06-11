% apply_dispersions.m
%   Apply uncertainties to input structure
%
% Aeroassist Simulation source code
%
% Inputs:
%   in: input structure
% 
% Output:
%   in: perturbed input structure
%   sigs: vector of standard deviations
%       1: m (kg)
%       2: cl
%       3: cd
%       4: aoa (rad)
%       5: vmag_atm (m/s)
%       6: hmag_atm (m)
%       7: fpa (rad)
%       8: lat (rad)
%       9: lon (rad)
%       10: az (rad)
% 
% Major Revision History:
%   *Created APR 2016, E. Roelke

function in = apply_dispersions(in,sigs)
%#codegen

%alter initial state vector to include random value for given mean and
%sigma value uncertainties
%% mass
sigma_m = in.v.mp.m_ini * sigs.m;
sigma_cd = in.v.aero.cd * sigs.cd;
if (in.v.gnc.g.p.dm_mode == uint8(4))
    scale = sqrt(cast(in.v.gnc.g.p_dej_n.n_jett + 1,'double') / 2);
    sigma_m = sigma_m / scale;
    sigma_cd = sigma_cd / scale;
end
in.v.mp.m_ini = in.v.mp.m_ini + normrnd(0, sigma_m);

%% Aerodynamics
% in.v.aero.cl = in.v.aero.cl + normrnd(0,s_cl);
in.v.aero.cd = in.v.aero.cd + normrnd(0, sigma_cd);
% in.v.aero.aoa = in.v.aero.aoa + normrnd(0,s_aoa);

%% Entry state
% v0 = norm(in.s.traj.v_pci_ini);
% % h0 = norm(in.s.traj.r_pci_ini)-in.p.r_e;
% h0 = in.s.traj.alt;     %initial altitude
% V_azi = asin(in.s.traj.v_pci_ini(3)/v0);
% fpa0 = asin(in.s.traj.v_pci_ini(1)/(v0*cos(V_azi)));

%Position uncertanties
% %x error (5.62 km 1 sigma)
% in.s.traj.r_pci_ini(1) = in.s.traj.r_pci_ini(1) + normrnd(0,5.62e3 * 3);
% %y error (32.1 km 1 sigma)
% in.s.traj.r_pci_ini(2) = in.s.traj.r_pci_ini(2) + normrnd(0,32.1e3 * 3);
% %z error (3.63 km 1 sigma)
% in.s.traj.r_pci_ini(3) = in.s.traj.r_pci_ini(3) + normrnd(0,3.63e3 * 3);

%velocity uncertainty (1-sigma)
in.s.traj.vel_pp_mag = in.s.traj.vel_pp_mag + normrnd(0,sigs.vmag_atm);

%FPA variation (1-sigma)
in.s.traj.gamma_pp = in.s.traj.gamma_pp + normrnd(0,sigs.efpa);

%Entry altitude variation (1-sigma)
in.s.traj.alt = in.s.traj.alt + normrnd(0,sigs.hmag_atm);

%Entry latitude and longitude (.12deg and .13deg (1 sigma), respectively)
% in.s.traj.lat = in.s.traj.lat + normrnd(0,sigs.lat0);
% in.s.traj.lon = in.s.traj.lon + normrnd(0,sigs.lon0);
% in.s.traj.az = in.s.traj.az + normrnd(0,sigs.az0);

% dej n uncertainties
if (in.v.gnc.g.p.atm.Kflag)
    switch (in.v.gnc.g.p.dm_mode)
        case 3 %pdg
            in.v.gnc.g.pd.sigs.m = sigs.m;
            in.v.gnc.g.pd.sigs.cd = sigs.cd;
            in.v.gnc.g.pd.sigs.aref = 0;
        case 4  %dej_n
            in.v.gnc.g.p_dej_n.sigs.m = sigs.m;
            in.v.gnc.g.p_dej_n.sigs.cd = sigs.cd;
            in.v.gnc.g.p_dej_n.sigs.aref = 0;
    end
end

% update inertial entry state
[in.s.traj.r_pci_ini, in.s.traj.v_pci_ini] = LLAVFA2RV_I( ...
    in.s.traj.lat, in.s.traj.lon, in.s.traj.alt, ...
    in.s.traj.gamma_pp, in.s.traj.az, in.s.traj.vel_pp_mag, ...
    in.p.omega, in.s.traj.t_ini, 0, in.p.r_e, in.p.r_p);
    
end