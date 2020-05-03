% update_t_jett.m
%   update jettison time estimate using hydra method
% 
% Entry systems Design Lab
%   University of Colorado Boulder
% 
% Written By: Evan Roelke, Aug 2019
% 
% Inputs:
%   tj_curr: current jettison time estimate (s)
%   dr1: apoapsis error estimate from tj_curr_i-1
%   dr2: apoapsis error estimate with tj_curr_i-1 + drErr
%   drErr: slight jettison time epsilon
%   dtj_lim: bound on jettison time update
%   p: dej_n guidance input struct
% 
function [tj_new, dtj] = update_t_jett(tj_curr, dr1, dr2, drErr, dtj_lim, tol_ap)
%#codegen

drdt = (dr2 - dr1)/drErr;
% d = f/f' in Newton Method (f = tgt apoapsis error)
dtj = dr1/drdt;

% put bounds on t_jett update (in case of hitting surface,
% hyperbolic orbit)
if dtj >= dtj_lim
    dtj = dtj_lim;
elseif dtj <= -dtj_lim
    dtj = -dtj_lim;
end %bounds

% if within tolerance and guidance tries to diverge, ignore the instability
if (abs(dtj) == dtj_lim && abs(dr1) <= tol_ap)
    dtj = 0;
    tj_new = tj_curr;
else
    tj_new = tj_curr - dtj;
end

end %update_t_jett