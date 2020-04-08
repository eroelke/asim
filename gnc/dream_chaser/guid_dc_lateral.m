% guid_dc_lateral.m
%   Dream Chaser entry guidance algorithm lateral logic
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   i - struct(1), multi, guidance input data structure
%   s - struct(1), multi, guidance state data structure
%   p - struct(1), multi, guidance parameter data structure
%
% Outputs:
%   s - struct(1), multi, guidance state data structure
%
% Major Revision History:
%   *Created 13 MAY, Z.R. Putnam, seperated from guidance.m

function [ s ] = guid_dc_lateral( i, s, p )
%#codegen

%% Compute lateral deadband                    
if s.vel_pp_mag > p.azi_limit_v1 
   s.azi_limit = p.azi_limit_1; % rad
elseif s.vel_pp_mag < p.azi_limit_v2
   s.azi_limit = p.azi_limit_2; % rad
else
   m = (p.azi_limit_2-p.azi_limit_3)/...
       (p.azi_limit_v2 - p.azi_limit_v1);
   s.azi_limit = p.azi_limit_3 - ...
       m*(p.azi_limit_v1 - s.vel_pp_mag);
end

%% Compute error
r_unit = [cos(s.lat)*cos(s.lon); ...
         cos(s.lat)*sin(s.lon); ...
         sin(s.lat)]; % nd
v_unit = i.vel_pp/s.vel_pp_mag; % nd
tgt_unit = [cos(p.latT)*cos(p.lonT);...
           cos(p.latT)*sin(p.lonT);...
           sin(p.latT)]; % nd
n_traj = cross(r_unit,v_unit); % nd
n_traj = n_traj/norm(n_traj);
n_tgt = cross(r_unit,tgt_unit); % nd
n_tgt = n_tgt/norm(n_tgt);
arg = median([-1, dot(n_traj,n_tgt),1]); % nd
if acos(dot(n_tgt,v_unit)) < pi/2
   s.azi_error = -acos(arg); % rad
else
   s.azi_error = acos(arg); % rad
end

%% Initialize if req'd
if i.init_flag
   s.cmd_bank_sign = int8(-sign(s.azi_error));
end

%% Check for bank reversal
if (double(s.cmd_bank_sign)*s.azi_error - s.azi_limit ) > double(0)
   % Initiate bank reversal 
   s.cmd_bank_sign = s.cmd_bank_sign*-1;
   s.bank_rev_count = s.bank_rev_count + uint8(1); % increment counter
end


end % guid_dc_lateral()