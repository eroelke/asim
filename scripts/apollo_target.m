%% Apollo targeting function
    % code xref:
    % guid_msl_prime.m	line 153

% inputs
%   gi: guidance input structure
%   gs: guidance state structure
%   gp: guidance parameter structure
%
% outputs:
%   gs: guidance state structure
%
%
%
%

function [ gs ] = apollo_target( gi, gs, gp )

% Compute guidance time
gs.time = gi.timetag - gs.time_init; % s

% TEST: has the relative velocity flag been set previously?
% if gs.V_rel_switch_reached == false
%   % NO, use inertial velocity
%   gs.V_work = gi.V_ocg_wrt_Eartho_icrf_ICRF; % m/s
%   
% else  
  % YES, use planet-relative velocity 
  gs.V_work = gi.V_ocg_wrt_Eartho_itrf_ICRF; % m/s
%   fprintf('switching to relative velocity @ t = %f, V = %f\n',gs.time,norm(gs.V_work,2))
  
% end


%% Compute vector magnitudes and quantities
% Velocity magnitude
gs.V_mag = norm(gs.V_work,2); % m/s
% Square of velocity magnitude, normalized by satellite velocity
gs.V_norm_sq = (gs.V_mag/gp.V_sat)^double(2); % nd
% Unit vector along current position vector
temp_v1 = gi.R_ocg_wrt_Eartho_ICRF / norm(gi.R_ocg_wrt_Eartho_ICRF,2); % nd
% Current altitude rate
gs.alt_rate = dot(gs.V_work, temp_v1); % m/s
% Vector normal to current trajectory plane
temp_v2 = cross(gs.V_work, temp_v1); % m/s
% Unit normal to current trajectory plane
gs.U_norm_traj = temp_v2 / norm(temp_v2,2); % nd
% Current acceleration magnitude, assumed to be aerodynamic only
gs.D_mag = norm(gi.Asens_onb_wrt_Eartho_icrf_OB,2); % m/s^2


%% Estimate target transfer angle
% TEST: has the relative velocity flag been set previously?
if gs.V_rel_switch_reached == false
  % NO, use inertial velocity

  % TEST: has the Final phase been reached?
%   if gs.phase < uint8(2)
%     % NO, Final phase not reached yet, 
%     % Use pre-Final phase method for target transfer angle calculation 
%     % (assumes a constant velocity flown to target)
%     gs.tgt_transfer_ang = gp.E_rotation_rate * ...
%       ( gp.tof_prm * gs.ang_to_tgt + gs.time ); % rad
% 
%   else   
    % YES, Final phase reached
    % Use inertial Final phase method for target transfer angle calculation 
    % (assumes current velocity flown all the way to target)
    gs.tgt_transfer_ang = gp.E_rotation_rate * ...
      ( gp.E_radius_eq * gs.ang_to_tgt / gs.V_mag + ...
      gs.time ); % rad

    % TEST: is velocity mag. low enough to switch to relative velocity?
    if (gs.V_mag - gp.V_rel_switch) <= double(0)
      % YES, switch to relative velocity
      gs.V_rel_switch_reached = true; % nd

    end
      % OTHERWISE, do not switch to relative velocity
      
%   end

else
  % YES, use relative velocity method for target transfer angle calculation (chase target)
  gs.tgt_transfer_ang = gp.E_rotation_rate * gs.time; % rad

end


%% Compute targeting angles and vectors
% Compute current target position unit vector
temp_s1 = cos( gs.tgt_transfer_ang ) - double(1); % nd
temp_s2 = sin( gs.tgt_transfer_ang ); % nd
gs.U_tgt = gs.U_tgt_init + gs.U_tgt_eq*temp_s1 + ...
	gs.U_tgt_east*temp_s2; % nd

% Ccompute lateral angle
gs.lat_ang = -dot( gs.U_tgt, gs.U_norm_traj ); % rad

% Compute longitudinal angle
gs.ang_to_tgt = acos( dot( gs.U_tgt, temp_v1 ) ); % rad

% Compute longitudinal range to target
gs.range_to_tgt = gs.ang_to_tgt * gp.E_radius_avg; % m

% determine next function call
switch gs.phase
    case 1 % Pre-bank Phase
        gs.next_step = uint8(1);
    case 2 % Final Phase
        gs.next_step = uint8(2);
    case 3 % Heading Alignment Phase
        gs.next_step = uint8(3);
    otherwise
        gs.next_step = uint8(5); % terminate guidance
end

end % apollo_target()

