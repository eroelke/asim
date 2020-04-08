% guid_dc_mach.m
%   Computes Mach number in Dream Chaser guidance
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   v_mag - double(1), m/s, wind-relative velocity magnitude
%   alt - double(1), m, altitude
%   p - struct(1), multi, guidance parameter data structure
%
% Outputs:
%   mach - double(1), nd, estimate of current Mach number
%
% Major Revision History:
%   *Created 5 APR 2012, Z.R.Putnam
%   *16 APR 2012, Z.R.Putnam, switched to guidance atmosphere model
%   *13 MAY 2012, Z.R.Putnam, removed "in" as an input

function mach = guid_dc_mach( v_mag, alt, p )
%#codegen

    T = lin_interp(p.atm_table(:,1),p.atm_table(:,4),alt); % K
    a = sqrt(p.atm_gamma*p.atm_R*T); % m/s
    mach = v_mag/a; % nd

end % guid_dc_mach()