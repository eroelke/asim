% get_dc_aero.m
% Determine trim AOA at input conditions using Golden Section algorithm
%   Ref: Vanderplaats, "Numerical Optimization Techniques for Engineering
%       Design: With Applications," 3rd edition, VR&D,Inc., 2001, p. 60.
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   act - struct(1), Dream Chaser actuator positions
%     DLE - double(1), deg, left elevon position
%     DRE - double(1), deg, right elevon position
%     DUL - double(1), deg, Upper left body flap position
%     DUR - double(1), deg, Upper right body flap position
%     DLL - double(1), deg, Lower left body flap position
%     DLR - double(1), deg, Lower right body flap position
%     DR - double(1), deg, Rudder position
%     GEAR - double(1), deg, gear position
%   beta - double(1), deg, sideslip angle
%   airspeed - double(1), ft/s, true airspeed
%   pqr - double(3), rad/s, body angular rate vector
%   agl - double(1), ft, altitude above surface
%   mach - double(1), nd, Mach number
%   qbar - double(1), psf, dynamic pressure
%   dc_aero - struct(1), md, Dream Chaser aerodynamics data
%   CG - double(3), ft, center of gravity location (CG)
%   CP - double(3), ft, center of pressure location (CP)
%   Sref - double(1), ft^2, aerodynamic reference area (Sref)
%   bref - double(1), ft, aerodynamic reference length (bref)
%   cbarref - double(1), ft, aerodynamic reference chord (cbarref)
%   
% Outputs:
%   alpha_trim - double(1), deg, trim angle-of-attack
%   My_trim - double(1), ft*lbf, pitching moment, M(2)
%
% Major Revision History:
%   *Created 15 MAR 2012, Z.R. Putnam

function [alpha_trim, My_trim] = get_dc_trim_aoa( act, beta, airspeed, ...
    pqr, agl, mach, qbar, dc_aero, dc_unc, CG, CP, Sref, bref, cbarref )
%#codegen

%% Parameters and settings
tau = 0.381966; % nd, golden ratio!
N = 20; % nd, number of function evaluations
ag_offset_high1 = 10; % deg, distance off nominal trim curve used to create initial bounds
ag_offset_high2 = 12; % deg, distance off nominal trim curve used to create initial bounds
ag_offset_low = -20; % deg, distance off nominal trim curve used to create initial bounds
mach_guess = [0,0.2,0.9,1.6,2.3,3.5,6,10,20.3,100]; % nd, Mach domain for initial guess curve
alpha_guess = [4.12,4.12,8.34,-1.30,2.21,23.40,34.14,33.36,37.65,37.65]; % deg, trim angle-of-attack curve for initial guess

if mach >= 5
    ag_offset_high = ag_offset_high1;
else
    ag_offset_high = ag_offset_high2;
end


%% Initialize

% Look up nominal angle-of-attack at input Mach number
ag = interp1(mach_guess,alpha_guess,mach);

% Lower bound
Xl = ag + ag_offset_low;
C = get_dc_aero_coeff( act, Xl, beta, airspeed, pqr, agl, mach, ...
    dc_aero, dc_unc, bref, cbarref );
[Fe,Me] = get_dc_aero_FM( C, qbar, Xl, beta, CG, CP, Sref, bref, cbarref ); %[lbf, ft*lbf]   
Fl = abs(Me(2));

% Upper bound
Xu = ag + ag_offset_high; % deg
C = get_dc_aero_coeff( act, Xu, beta, airspeed, pqr, agl, mach, ...
    dc_aero, dc_unc, bref, cbarref ); % nd
[Fe,Me] = get_dc_aero_FM( C, qbar, Xu, beta, CG, CP, Sref, bref, cbarref ); %[lbf, ft*lbf]   
Fu = abs(Me(2)); % ft*lbf

% Point 1
X1 = (1-tau)*Xl + tau*Xu; % deg
C = get_dc_aero_coeff( act, X1, beta, airspeed, pqr, agl, mach, ...
    dc_aero, dc_unc, bref, cbarref ); % nd
[Fe,Me] = get_dc_aero_FM( C, qbar, X1, beta, CG, CP, Sref, bref, cbarref ); %[lbf, ft*lbf]   
F1 = abs(Me(2)); % ft*lbf

% Point 2
X2 = tau*Xl + (1-tau)*Xu; % deg
C = get_dc_aero_coeff( act, X2, beta, airspeed, pqr, agl, mach, ...
    dc_aero, dc_unc, bref, cbarref ); % nd
[Fe,Me] = get_dc_aero_FM( C, qbar, X2, beta, CG, CP, Sref, bref, cbarref ); %[lbf, ft*lbf]   
F2 = abs(Me(2)); % ft*lbf

% Initialize best guess
bX = X1;
bF = F1;


%% Converge solution
kk = 3;
kk = kk + 1;
while(kk < N)

    if F1 > F2
        Xl = X1;
        Fl = F1;
        X1 = X2;
        F1 = F2;
        X2 = tau*Xl + (1-tau)*Xu;
        C = get_dc_aero_coeff( act, X2, beta, airspeed, pqr, agl, mach, ...
            dc_aero, dc_unc, bref, cbarref );
        [Fe,Me] = get_dc_aero_FM( C, qbar, X2, beta, CG, CP, Sref, bref, cbarref ); %[lbf, ft*lbf]   
        F2 = abs(Me(2));
        bX = X1; % current best X
        bF = F1; % current best F
    else
        Xu = X2;
        Fu = F2;
        X2 = X1;
        F2 = F1;
        X1 = (1-tau)*Xl + tau*Xu;
        C = get_dc_aero_coeff( act, X1, beta, airspeed, pqr, agl, mach, ...
            dc_aero, dc_unc, bref, cbarref );
        [Fe,Me] = get_dc_aero_FM( C, qbar, X1, beta, CG, CP, Sref, bref, cbarref ); %[lbf, ft*lbf]   
        F1 = abs(Me(2));
        bX = X2; % current best X
        bF = F2; % current best F
    end

    kk = kk + 1;

end


%% Assign output
alpha_trim = bX;
My_trim = bF;

end % get_dc_trim_aoa
