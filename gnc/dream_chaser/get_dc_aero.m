% get_dc_aero.m
%   Determine Dream Chaser aerodynamic coefficients, force, and moment for
%   several actuator options
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   aoa - double(1), rad, angle-of-attack
%   ssa - double(1), rad, sideslip angle
%   trim_flag - uint8(1), nd, trim angle-of-attack calculation flag
%   alt - double(1), m, altitude above surface
%   v_wind_mag - double(1), m/s, true airspeed
%   mach - double(1), nd, Mach number
%   dynp - double(1), Pa, dynamic pressure
%   cg - double(3), m, current c.g. position
%   in - struct(1), md, input data structure
%     from in.c: unit conversions
%     from in.v.aero:
%       dc_mode - uint8(1), nd, Dream Chaser aerodynamics mode
%       dc - struct(1), md, Dream Chaser aerodynamics data
%       r_cp - double(3), m, center of pressure location (CP)
%       area_ref - double(1), m^2, aerodynamic reference area (Sref)
%       span_ref - double(1), m, aerodynamic reference length (bref)
%       chord_ref - double(1), m, aerodynamic reference chord (cbarref)
%   
% Outputs:
%   C - double(6), nd, aerodynamic coefficients
%      C(1) - Cd - drag coefficient
%      C(2) - Cybeta
%      C(3) - Cl - lift coefficient
%      C(4) - Clbeta
%      C(5) - Cm_ref
%      C(6) - Cnbeta_ref
%   F - double(3), N, aerodynamic force
%   M - double(3), N*m, aerodynamic moment
%   aoa - double(3), 
%   act - struct(1), Dream Chaser actuator positions
%     DLE - double(1), deg, left elevon position
%     DRE - double(1), deg, right elevon position
%     DUL - double(1), deg, Upper left body flap position
%     DUR - double(1), deg, Upper right body flap position
%     DLL - double(1), deg, Lower left body flap position
%     DLR - double(1), deg, Lower right body flap position
%     DR - double(1), deg, Rudder position
%     GEAR - double(1), deg, gear position
%
% Major Revision History:
%   *Created 22 FEB 2012, Z.R. Putnam
%   *16 MAR 2012, Z.R. Putnam, updated to include addition options,
%       integrated into simulation

function [C,F,M,aoa,act] = get_dc_aero( aoa, ssa, trim_flag, ...
   alt, v_wind_mag, mach, dynp, cg, in )
%#codegen

%% Initialize
% Initialize output
C = nan(6,1); % nd, aerodynamic coefficients
F = nan(3,1); % lbf, aerodynamic force
M = nan(3,1); % ft*lbf, aerodynamic moment
act = in.v.aero.dc_act; % deg, DC steady-state actuator deflections
% Convert units to English
alpha = aoa*in.c.r2d; % deg, convert from rad
beta = ssa*in.c.r2d; % deg, convert from rad
agl = alt*in.c.m2ft; % ft, convert from m
airspeed = v_wind_mag*in.c.m2ft; % ft/s, convert from m/s
qbar = dynp*in.c.Pa2psf; % psf, convert from Pa
CG = cg*in.c.m2ft; % ft, convert from m
CP = in.v.aero.r_cp*in.c.m2ft; % ft, convert from m
Sref = in.v.aero.area_ref*in.c.m2ft*in.c.m2ft; % ft^2, convert from m^2
bref = in.v.aero.span_ref*in.c.m2ft; % ft, convert from m
cbarref = in.v.aero.chord_ref*in.c.m2ft; % ft, convert from m     
pqr = [0;0;0]; % rad/s, body angular rates

%% Compute aerodynamics
switch in.v.aero.dc_mode
         
    case 1
    %% No actuator deflection - forced angle-of-attack
    % Compute coefficients, force, moment at given aoa, ssa, state
    
        C = get_dc_aero_coeff( act, alpha, beta, airspeed, pqr, agl, ...
            mach, in.v.aero.dc, in.v.aero.dc_unc, bref, cbarref );
        [Fe,Me] = get_dc_aero_FM( C, qbar, alpha, beta, ...
            CG, CP, Sref, bref, cbarref ); %[lbf, ft*lbf]
        
        
    case 2
    %% No actuator deflections - trim angle-of-attack
    % Compute coefficients, force, moment, trim aoa at given ssa, state
        
        if trim_flag
            % Determine trim angle-of-attack
            [alpha,My] = get_dc_trim_aoa( act, beta, airspeed, pqr, agl, ...
                mach, qbar, in.v.aero.dc, in.v.aero.dc_unc, CG, CP, Sref, ...
                bref, cbarref ); % [deg,ft*lbf]
        end
        
        % Get aerodynamic coefficients and forces
        C = get_dc_aero_coeff( act, alpha, beta, airspeed, pqr, agl, ...
            mach, in.v.aero.dc, in.v.aero.dc_unc, bref, cbarref ); % nd
        [Fe,Me] = get_dc_aero_FM( C, qbar, alpha, beta, ...
            CG, CP, Sref, bref, cbarref ); % [lbf, ft*lbf]
        
        
    case 3
    %% Bank-rate command - trim angle-of-attack
    % Compute coefficients, force, moment, trim aoa, elevon deflection at
    % given ssa, state
    
        C = zeros(6,1);
        Fe = zeros(3,1); % lbf
        Me = zeros(3,1); % ft*lbf
        
        
    otherwise
    %% Incorrect mode option, set output to zero
        
        C = zeros(6,1);
        Fe = zeros(3,1); % lbf
        Me = zeros(3,1); % ft*lbf

        
end % switch dc_act_mode


%% Convert outputs to metric units
F = Fe/in.c.N2lbf; % lbf -> N
M = Me/in.c.N2lbf/in.c.m2ft; % ft*lbf -> N*m
aoa = alpha/in.c.r2d; % deg -> rad


end % get_dc_aero
 
