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
in.v.mp.m_ini = in.v.mp.m_ini + normrnd(0,sigs.m);

%% Aerodynamics
% in.v.aero.cl = in.v.aero.cl + normrnd(0,s_cl);
in.v.aero.cd = in.v.aero.cd + normrnd(0,sigs.cd);
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

%% Atmosphere
%multiplicative factor on atmospheric properties
% in.p.atm.table(:,1) = in.p.atm.table(:,1).*disp_atmos(1); %density
% in.p.atm.table(:,2) = in.p.atm.table(:,2).*disp_atmos(2); %temperature
% in.p.atm.table(:,3) = in.p.atm.table(:,3).*disp_atmos(3); %pressure
% in.p.atm.table(:,4) = in.p.atm.table(:,4).*disp_atmos(4); %easterly winds
% in.p.atm.table(:,5) = in.p.atm.table(:,5).*disp_atmos(5); %northern winds
% in.p.atm.table(:,6) = in.p.atm.table(:,6).*disp_atmos(6); %vertical winds

%% Deceleretor uncertainties
%chute deploy altitude

%chute deploy velocity

%chute aerodynamics (cd +- 0.01 1 sigma) - random guess
% in.v.decel.cd = in.v.decel.cd + normrnd(0,0.01);

% jettisoned mass
% in.v.gnc.g.p_sej_a.mass_jettison = in.v.gnc.g.p_sej_a.mass_jettison + normrnd(0,sigs(1) );


end