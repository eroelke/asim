% guid_dc_pred.m 
%   Predictor for Dream Chaser entry guidance algorithm
%
% Aeroassist Simulation source code
%
% Developed by:
%   Space Systems Design Lab
%   Daniel Guggenheim School of Aerospace Engineering
%   Georgia Institute of Technology
%
% Inputs:
%   y0 - double(4), multi, initial state vector for propogation
%   t0 - double(1), s, initial time
%   clsf - double(1), nd, lift coefficient scale factor
%   cdsf - double(1), nd, drag coefficient scale factor
%   p - struct(1), multi, guidance parameter structure
%   bank - double(1), rad, current bank command
%   r_vec - double(3), m, initial position vector
%   v_vec - double(3), m/s, initial velocity vector
%   t_max - double(1), s, maximum allowable integration time
%   
% Outputs:
%   t - double(1), s, final time
%   y - double(4), multi, final state
%   flag - int8(1), nd, termination flag
%   maxHeating - double(1), W/m^2, peak heat rate over integration
%   maxGloading - double(1), nd, peak g-load over integration
%
% Major Revision History:
%   *Created 3 FEB 2012, Z.R. Putnam
%   *16 JUL 2012, Z.R. Putnam, updated comments and arguments

function [t,y,flag,maxHeating,maxGloading] = guid_dc_pred(y0, t0, clsf, ...
    cdsf, p, bank, r_vec, v_vec, t_max)
%#codegen

%% Initialize
y = y0; % multi, initial state
t = t0; % s, initial time
area_ref = p.area_ref; % m^2, set reference area to initial value
flag = 0; % nd, termination flag
maxHeating = 0; % W/cm^2, maximum heat rate observed during prediction
qDotPrior = 0; % W/cm^2, heat rate of prior time step
maxGloading = 0; % Earth-g, maximum g-loading observed during propagation
nPrior = 0; % Earth-g, g-loading of prior time step
checkHeating = true; % int8, flag to determine if maximum heat rate should be checked
observedHeatingIncrease = false; % int8, flag to determine if heat rate increase has been observed
checkGLoading = true; % int8, flag to determine if maximum g-loading should be checked
observedGLoadingIncrease = false; % int8, flag to determine if g-loading increase has been observed

%% Integration (RK4)
while ~flag
   
    %% Integrate one time step
    [k1,L,D] = eom_planar(y, area_ref, clsf, cdsf,p, bank, r_vec, v_vec);
    [k2,~,~] = eom_planar(y+0.5*p.h*k1, area_ref, clsf, cdsf, p, bank, r_vec, v_vec);
    [k3,~,~] = eom_planar(y+0.5*p.h*k2, area_ref, clsf, cdsf,p, bank, r_vec, v_vec);
    [k4,~,~] = eom_planar(y+p.h*k3, area_ref, clsf, cdsf, p, bank, r_vec, v_vec);
    y = y + (p.h/6)*(k1 + 2*k2 + 2*k3 + k4);
    t = t + p.h; % increment time
    
    %% Update dependent states
    alt = y(1) - p.p_r;
    v = y(3);

    %% Compute maximum constraint values
    
    % Heating
    vel_pp_mag = y(3);
    % Get density
    rho = lin_interp(p.atm_table(:,1), p.atm_table(:,2), alt);
    % Compute heat rate
    qDot = p.k*sqrt(rho/p.nose_radius)*vel_pp_mag^3;    
    
    % Record maximum heat rate if recording functionality activated through
    % checkHeating variable.
    if (qDot > maxHeating) && checkHeating
      maxHeating = qDot;
    end
    
    % Determine if heating increase has been observed.
    if (qDot > qDotPrior) && ~observedHeatingIncrease
      observedHeatingIncrease = true;
    end
    
    % If heating increase observed and heat rate is decreasing, stop recording
    % maximum heat rate.
    if (qDot < qDotPrior) && observedHeatingIncrease
      checkHeating = false;
    end
    qDotPrior = qDot;
    
    % G-loading
    n = sqrt(L^2+D^2)/p.mass/p.surface_g;
    
    % Record maximum g-loading if recording functionality activated through
    % checkGLoading variable.
    if (n > maxGloading) && checkGLoading
      maxGloading = n;
    end
    
    % Determine if g-loading increase has been observed.
    if (n > nPrior) && ~observedGLoadingIncrease
      observedGLoadingIncrease = true;
    end
    
    % If g-loading increase observed and g-loading is decreasing, stop recording
    % maximum g-loading.
    if (n < nPrior) && observedGLoadingIncrease
      checkGLoading = false;
    end
    nPrior = n;
    
    %% Check termination conditions
    if alt > p.alt_max
      flag = 1; % Maximum altitude exceeded
      break;
    elseif alt < p.alt_min
      flag = 2; % Minimum altitude exceeded
      break;
    elseif t > t_max
      flag = 3; % Maximum time exceeded
      break;
    elseif v < p.vMin
      flag = 4; % Minimum velocity exceeded
      break;
    end
end

end % guid_dc_pred

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ydot,L,D] = eom_planar(y, area_ref, clsf, cdsf, p, bank, r_vec, v_vec)
% Predictor function. Equations of motion in plane.

%% Assign state vectors
r = y(1); % m, radius from center of Earth
%theta = y(2); % rad, downrange angle
v = y(3); % m/s, relative velocity magnitude
gam = y(4); % rad, relative flight path angle

%% Compute atmospheric properties
alt = r - p.p_r; % m, geocentric altitude
dens = lin_interp(p.atm_table(:,1), p.atm_table(:,2),alt); % kg/m^3, atmosphere model density
T = lin_interp(p.atm_table(:,1),p.atm_table(:,4),alt); % K, atmospheric temperature
a = sqrt(p.atm_gamma*p.atm_R*T); % m/s, speed of sound
mach = v/a;

%% Compute aerodynamic coefficients as a function of nominal angle of attack profile
[cl, cd, lod] = guid_dc_get_est_aero( clsf, cdsf, p.aero_model, mach);

%% Aerodynamic forces
dynp = 0.5*dens*v^2; % N/m^2, dynamic pressure
D = cd*dynp*area_ref; % N, drag force
L = cl*dynp*area_ref; % N, vertical lift force

%% Coriolis force
% This force is usually not included within the planar equations of motion
% since it cannot be calculated directly from the planar states alone. By
% providing the current planet-relative velocity vector direction, the Coriolis
% force is approximated assuming the direction of the velocity vector does not
% change during the propagation. Expansion of the predictor to six states would
% allow improved calculations during the propagation. However, this
% approximation is currently sufficient. Without this approximation, the
% downrange predictions would not be accurate.
F3D = p.mass*2*cross(p.omega,v_vec/norm(v_vec)*v); % N, approximate Coriolis force magnitude
F = -dot(F3D,r_vec/norm(r_vec))/cos(gam); % N, component of approximate Coriolis force in radial direction

%% State derivatives
rDot = v*sin(gam); % m/s, d(radius)/d(time)
thetaDot = v*cos(gam)/r; % rad/s, d(downrange angle)/d(time)
vDot = -D/p.mass - p.mu*sin(gam)/r^2; % m/s^2, d(relative velocity magnitude)/d(time)
gamDot = (L*cos(bank)+F)/(p.mass*v) + (v/r - p.mu/(v*r^2))*cos(gam); % rad/s, d(relative flight path angle)/d(time)


%% Output
ydot = [rDot; thetaDot; vDot; gamDot]; % state derivative vector

end % eom_planar

