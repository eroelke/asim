% guid_gt.m 
%   Numeric prediction algorithm to determine when to begin constant thrust
%   gravity turn.
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   i - struct(1), multi, guidance input structure
%   s - struct(1), multi, guidance state structure
%   p - struct(1), multi, guidance parameter structure
%   init_flag - logical(1), nd, first-pass initialization flag
%   
% Outputs:
%   s - struct(1), multi, guidance state structure
%
% Major Revision History:
%   *Created 11 NOV 2011, Z.R. Putnam
%   *1 AUG 2012, Z.R.Putnam, updated, removed rate limiter

function [s] = guid_gt(i,s,p,init_flag)
%#codegen

%% First-pass initialization
if init_flag
    s.ignition = false;
end

%% Check for ignition time
V_mag = norm(i.V_pf); % m/s, current planet-fixed velocity magnitude

if (~s.ignition) && (V_mag < p.V_mag_min)

    %% Predict trajectory to terminal velocity based on current state
    v0 = V_mag; % m/s, velocity magnitude
    h0 = norm(i.R) - p.r_e; % m, altitude

    uR = i.R/norm(i.R); % nd, unit vector along position vector
    uV = i.V_pf/norm(i.V_pf); % nd, unit vector along velocity
    psi0 = acos(dot(-uR,uV)); % rad, psi angle between gravity and planet-fixed velocity

    y0 = [v0; psi0; h0; p.m_ini]; % multi, initial state vector

    
    [y,flag] = predictor(y0,p); %  multi, propagate trajectory to velocity target
    %% If predicted terminal altitude at or below target, set ignition flag
    switch flag
        case 1 % target velocity reached
            if y(3) < p.tgt_alt
                s.ignition = true;
            end
        case 2 % negative mass
            % do nothing!
        case 3 % time out
            % do nothing!
        otherwise
            % do nothing!
    end % switch flag

end

end % guid_gt



function [y,flag] = predictor(y0, p)
%% Propagate equations of motion to velocity target

%% Initialize
y = y0; % multi, current state
V_error = 0; % m/s, current velocity error
V_error_prev = 0; % m/s, previous velocity error
t = 0; % s, initial time
flag = 0; % nd, termination flag

%% Integration (RK4)
while ~flag
    
    %% Integrate one time step
    k1 = eom(y, p);
    k2 = eom(y+0.5*p.h*k1, p);
    k3 = eom(y+0.5*p.h*k2, p);
    k4 = eom(y+p.h*k3, p);
    y = y + (p.h/6)*(k1 + 2*k2 + 2*k3 + k4);
    t = t + p.h; % increment time
       
    %% Check termination conditions
    V_error = y(1)-p.tgt_V_mag;
    if V_error*V_error_prev < 0
        flag = 1; % zero crossed
        break;
    elseif y(4) <= 0
        flag = 2; % Negative mass
        break;
    elseif t >= p.t_max
        flag = 3; % Maximum time exceeded
        p.t_max
        break;
    end
    
    V_error_prev = V_error;
end
    
end % predictor



function [ydot] = eom(y, p)
%% Equations of motion for predictor

% Assign state
v = y(1); % velocity
psi = y(2); % psi angle
h = y(3); % altitude
m = y(4); % mass

vdot = p.g*cos(psi) - p.thrust/m; % m/s^2
psidot = -p.g/v * sin(psi); % rad/s
hdot = -v*cos(psi); % m/s
mdot = -p.thrust/(p.g0*p.isp); % kg/s

% Compile ydot
ydot = [vdot; psidot; hdot; mdot];

end % eom

