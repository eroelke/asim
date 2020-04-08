% navigation_init.m
%   Initialize navigation model
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   in - struct(1), multi, input data structure
%
% Outputs:
%   nav - struct(1), multi, navigation data structure
%
% Major Revision History:
%   *Created 9 AUG 2012, Z.R.Putnam

function [ nav, nav_rnd ] = navigation_init( in )

%% Construct navigation data structure
nav = define_nav;

%% Copy navigation parameters
nav.p = in.v.gnc.n.p;
if nav.p.mode ~= 2
    nav.p.omega = in.p.omega;
    nav.p.r_e = in.p.r_e;
    nav.p.r_p = in.p.r_p;
end

%% Initialize random numbers for navigation error model
t_length = round((in.s.traj.t_max-in.s.traj.t_ini)*in.s.traj.rate)+1;
if nav.p.mode == 2
%     rng(nav.p.seed,'combRecursive'); % use input seed
%     rng('shuffle')
    nav_rnd = randn(9,t_length); % create repeatable set of random numbers on U[-1,1] for navigation model
    nav.s.x_ercv = (1 - exp(-2*nav.s.dt/nav.p.tau))*nav_rnd(:,1);
    nav.s.rva_error = nav.p.P_SS*nav.s.x_ercv;
else
    nav_rnd = zeros(9,t_length+1);
    nav.s.x_ercv = zeros(9,1);
    nav.s.rva_error = zeros(9,1);
end

%% Set states
nav.s.r_pci = in.s.traj.r_pci_ini + nav.s.rva_error(1:3);
nav.s.v_inrtl_pci = in.s.traj.v_pci_ini + nav.s.rva_error(4:6);
nav.s.a_sens_pci = zeros(3,1) + nav.s.rva_error(7:9);
% Set derived states to zero
nav.s.r_pcpf = zeros(3,1);
nav.s.v_pf_pci = zeros(3,1);
nav.s.v_pf_pcpf = zeros(3,1);
nav.s.a_sens_pcpf = zeros(3,1);
nav.s.alt = 0;
% Time and rate parameters
nav.s.t_ini = in.s.traj.t_ini;
nav.s.dt = 1/nav.p.rate; % seconds
nav.s.t = -nav.s.dt; % start one time step back
if nav.p.rate == 0
    nav.s.sn_ratio = 1; % execute every sim time step
else
    nav.s.sn_ratio = round(in.s.traj.rate/nav.p.rate); % sim/navigation rate ratio, nd
end

end % naviation_init()
