% control_init.m 
%   Initializes control data structure
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   in - data structure, misc, input data structure
%   
% Outputs:
%   ctrl - data structure, misc, control data structure with parameter
%       values set
%
% Major Revision History:
%   *Created 17 NOV 2011, Z.R. Putnam

function [ctrl] = control_init( in )
%#codegen

%% Create data structure
ctrl = define_ctrl;

%% Initialize parameters
ctrl.p = in.v.gnc.c.p;

%% Initialize states
% Time state
ctrl.s.t_ini = in.s.traj.t_ini; % s
ctrl.s.t_prev = in.s.traj.t_ini; % s
ctrl.s.t = in.s.traj.t_ini; % s
ctrl.s.dt = 0; % s
ctrl.s.bank_accel = 0; % rad/s^2
ctrl.s.bank_rate = in.v.gnc.g.p.const_bank_rate; % rad/s
ctrl.s.bank = in.v.gnc.g.p.const_bank; % rad
ctrl.s.area_ref = in.v.aero.area_ref; % m^2
ctrl.s.sc_ratio = round(in.s.traj.rate/ctrl.p.rate); % sim/control rate ratio, nd
ctrl.s.cd = in.v.aero.cd;

ctrl.s.jettison = false;
ctrl.s.ignition = false;
ctrl.s.deploy_decel = false;

end % control_init