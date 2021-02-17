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
nav.p.atm.omega = in.p.omega;
nav.p.atm.r_e = in.p.r_e;
nav.p.atm.r_p = in.p.r_p;

%% Initialize random numbers for navigation error model
nav.s.rva_ind = 1;
t_length = round((in.s.traj.t_max-in.s.traj.t_ini)*in.s.traj.rate)+1;
nav_rnd = randn(9,t_length); % preallocate randomly-sampled standard gaussian

switch nav.p.mode
    case 2 % ecrv
    %     rng(nav.p.seed,'combRecursive'); % use input seed
    %     rng('shuffle')
        nav.s.x_ercv = (1 - exp(-2*nav.s.dt/nav.p.ecrv.tau))*nav_rnd(:,1);
        nav.s.rva_error = nav.p.ecrv.P_SS*nav.s.x_ercv;
    case 5 % exponentially-decaying ecrv
        P = nav.p.n_exp.P0;

        nav.s.ercv = nav.p.n_exp.ercv0;%(1 - exp(-nav.s.dt/nav.p.ecrv.tau))*nav_rnd(:,1);    %n_noise
        nav.s.x_ercv = nav.s.ercv;
        nav.s.rva_error = real(sqrtm(P))*nav.s.x_ercv;  
    case 6 % predetermined exp-decay ecrv
        nav.s.rva_error = nav.p.rva_errs(:,1);
        nav.s.rva_ind = 2;
    otherwise
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
