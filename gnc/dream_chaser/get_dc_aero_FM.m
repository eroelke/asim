% get_dc_aero_FM.m 
%   Determine Dream Chaser aerodynamic forces and moments; based on
%   Simulink model provide by Ernie Lagimoniere (SNC)
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%  
% Based on:
%   Sierra-Nevada Corp. Dream Chaser aero_5p2_standalone.mdl
%
% Inputs:
%   C_total - double(6), nd, aerodynamic coefficients
%      C_total(1) - Cd - drag coefficient
%      C_total(2) - Cybeta
%      C_total(3) - Cl - lift coefficient
%      C_total(4) - Clbeta
%      C_total(5) - Cm_ref
%      C_total(6) - Cnbeta_ref
%   qbar - double(1), psf, dynamic pressure
%   alpha - double(1), deg, angle-of-attack
%   beta - double(1), deg, sideslip angle
%   CG - double(3), ft, center of gravity location
%   CP - double(3), ft, center of pressure location
%   Sref - double(1), ft^2, aerodynamic reference area
%   bref - double(1), ft, aerodynamic reference length
%   cbarref - double(1), ft, aerodynamic reference chord
%   
% Outputs:
%   F - double(3), lbf, aerodynamic force, body frame, body coordinates
%   M - double(3), ft*lbf, aerodynamic moment about CG, body frame, body
%      coordinates
%
% Major Revision History:
%   *Created 21 FEB 2012, Z.R. Putnam

function [F,M] = get_dc_aero_FM( C_total, qbar, alpha, beta, CG, CP, ...
    Sref, bref, cbarref )
%#codegen

%% Force

CF = [-C_total(1); C_total(2); -C_total(3)];

% Coeffecient transformation

% Incidence, Sideslip, & Airspeed
aoa = alpha*(pi/180); % rad, convert angle-of-attack to radians
ssa = beta*(pi/180); % rad, assume zero sideslip angle

% Direction Cosine Matrix Body to Wind
sa = sin(aoa);
sb = 0;
ca = cos(aoa);
cb = 1;
DCM = transpose([ca*cb,sb,sa*cb; -sb*ca,cb,-sa*sb; -sa,0,ca]);
Cb = DCM*CF;

F = qbar*Sref*Cb;


%% Moment

CM = C_total(4:6);
CGrCP = CG - CP;

% Subsystem
Out1 = [(CGrCP(2)*Cb(3) - CGrCP(3)*Cb(2))/bref; ...
        (CGrCP(3)*Cb(1) - CGrCP(1)*Cb(3))/cbarref; ...
        (CGrCP(1)*Cb(2) - CGrCP(2)*Cb(1))/bref ];

CM = CM + Out1;
qSL = qbar*Sref*[bref; cbarref; bref];
M = qSL .* CM;


end % get_dc_aero_FM