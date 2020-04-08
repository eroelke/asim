% guid_dc_main.m
%   Dream Chaser entry guidance algorithm
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
%   init_flag - logical(1), nd, initialization flag (true for first call to
%       guidance)
%
% Outputs:
%   s - struct(1), multi, guidance state data structure
%
% Major Revision History:
%   *Created 13 MAY 2012, Z.R. Putnam, seperated from guidance.m

function [ s ] = guid_dc_main( i, s, p )
%#codegen

%% First-pass initialization
if i.init_flag
    % Initialize aerodynamic estimator
    s.clsf = 0; % nd
    s.cdsf = 0; % nd

    % Set number of bank reverals to zero
    s.bank_rev_count = uint8(0); % nd
    
    % Set heat rate controller states
    s.heatRateDelayTrigger = int8(0); % nd
    s.heatRateDelayTime = 5; % s
    s.heatRateTriggerTime = p.t_max_ini; % s
    
    % Other initialization
    s.t_max = p.t_max_ini; % s
    
    % modes
    s.pcMode = int8(1); % nd
    s.egMode = int8(3); % nd
    
end % init_flag


%% Initialize at current time step
cmd_bank_prior = s.cmd_bank;

% Compute derived quantities
s.pos_pp_mag = norm(i.pos_pp,2); % position vector magnitude, m
s.vel_pp_mag = norm(i.vel_pp,2); % planet-relative velocity vector magnitude, m/s
u_pos_pp = i.pos_pp/s.pos_pp_mag; % unit vector along pos_pp, nd
u_vel_pp = i.vel_pp/s.vel_pp_mag; % unit vector along vel_pp, nd
arg = median([-1,1,dot(u_pos_pp, u_vel_pp)]); % divide-by-zero protection
s.gamma_pp = pi/2 - acos(arg); % planet-relative flight-path angle, rad

% Compute latitude and longitude (geocentric)
s.lat = asin(i.pos_pp(3)/s.pos_pp_mag);
s.lon = asin(i.pos_pp(2)/sqrt(i.pos_pp(1)^2+i.pos_pp(2)^2));
% Quadrant check on longitude
if i.pos_pp(1) < 0
    if i.pos_pp(2) > 0
        s.lon = pi - s.lon;
    elseif i.pos_pp(2) < 0
        s.lon = -pi - s.lon;
    end
end

if norm(i.sens_accel_ii) > 0.05
    % Compute Mach number
    s.mach = guid_dc_mach( s.vel_pp_mag, i.alt, p );
    % Compute planet-relative velocity vector in inertial coordinate system               
    v_rel_ii = i.vel_ii - cross(p.omega,i.pos_ii);

    [s.cd_est, s.cl_est, s.cdsf, s.clsf] = guid_dc_aero_est( ...
        v_rel_ii, i.sens_accel_ii, s.mach, p.aero_model, ...
        i.alt, s.clsf, s.cdsf, p );

    % Update velocity gate: ensures range targeting doesn't happen too 
    %   early; limits range targeting to velocities below vEndQdot
    [cl,cd,lod] = guid_dc_get_est_aero( s.clsf, s.cdsf, p.aero_model, 25); % Get aerodynamic properties at Mach 25
    if cd > p.vEndQDot_cd_lim
        s.vEndQDot = p.vEndQDot_high_cd; % m/s
    else
        s.vEndQDot = p.vEndQDot_low_cd; % m/s
    end

end


%% Main guidance modes
switch s.pcMode

    case 1 % Pullout
      
      % Fly full lift up through pullout
      s.cmd_bank = 0; % rad

      % When flight path angle levels off, switch to downrange targeting mode
      if s.gamma_pp > p.caGamma % consider alternate parameter such as g-loading

        % Change to energy depletion mode
        s.pcMode = int8(2);

        % Initialize initial location connected to target location for
        % ranging calculations using great circle.
        s.lat0 = s.lat;
        s.lon0 = s.lon;

      end

    case 2 % Downrange targeting mode

           
      % Run range calculation
      s.egMode = int8(3); % Bisection search to calculate bank to converge downrange to target.
      s = guid_dc_npc( i, s, p );
      cmd_bank_range = s.cmd_bank; % Bank command to converge to downrange target, rad

      % Run heat rate calculation
      s.egMode = int8(4); % Fly to heat rate constraint
      s = guid_dc_npc( i, s, p );
      cmd_bank_heating = s.cmd_bank; % Bank command to fly to heat rate constraint, rad

      % Run g-loading calculation
      s.egMode = int8(5); % Fly to g-loading constraint
      s = guid_dc_npc( i, s, p );
      cmd_bank_gloading = s.cmd_bank; % Bank command to fly to g-loading constraint, rad

      % Determine bank that satisfies constraints
      s.cmd_bank = min(cmd_bank_heating, cmd_bank_gloading);

      % Limit bank angle extremes early in entry to minimize large changes in
      % control that would otherwise be commanded during low dynamic pressure
      % regimes.
      if s.vel_pp_mag > 22000*0.3048
        if s.cmd_bank > 150*pi/180
          s.cmd_bank = 150*pi/180;
        elseif s.cmd_bank < 30*pi/180
          s.cmd_bank = 30*pi/180;
        end
      end

      % If bank angle to converge range is less than the bank required to
      % satisfy in-flight constraints, perform range targeting. Also ensure
      % velocity is below vEndQDot so range targeting does not occur too early.
      if cmd_bank_range + p.rangeBankBias < s.cmd_bank && ...
         s.vel_pp_mag < s.vEndQDot    

        % To ensure enough energy is depleted for high L/D cases, only perform
        % range targeting after a timed delay.
        
        % Determine if heat rate delay trigger activated
        if ~s.heatRateDelayTrigger

          % Initialize heat rate delay values
          s.heatRateDelayTrigger = int8(1);
          s.heatRateTriggerTime = i.t;

          % Compute the timed delay
          [cl,cd,lod] = guid_dc_get_est_aero( s.clsf, s.cdsf, p.aero_model, 25);
          if lod < 1.05
              s.heatRateDelayTime = 5; % s
          else
              s.heatRateDelayTime = (40/0.15)*(lod-1.05) + 5; % s
          end

        elseif i.t > s.heatRateTriggerTime + s.heatRateDelayTime

          % Timer delay has finished. Perform downrange targeting.
          s.cmd_bank = cmd_bank_range;

        end
        
      end

      % Activate heading alignment. Limit maximum bank angle.
      if s.vel_pp_mag < p.haVel
        s.pcMode = int8(3);
      end

    case 3
      % Heading alignment. Limit maximum bank angle.
      if abs(s.cmd_bank) > p.haMaxBank
        s.cmd_bank = sign(s.cmd_bank)*p.haMaxBank;
      end

  otherwise

    s.cmd_bank = 0; % rad

end


%% Lateral guidance
s = guid_dc_lateral( i, s, p );


%% Assign and filter bank command
s.cmd_bank = abs(s.cmd_bank)*double(s.cmd_bank_sign); % rad

s.cmd_bank = (1-p.bankFilterConst)*cmd_bank_prior + ...
    (p.bankFilterConst)*s.cmd_bank; % rad


end % guid_dc_main()
